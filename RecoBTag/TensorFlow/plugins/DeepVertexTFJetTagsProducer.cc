
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/DeepFlavourTagInfo.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "RecoBTag/TensorFlow/interface/tensor_fillers.h"


// Declaration of the data structure that is hold by the edm::GlobalCache.
// In TensorFlow, the computational graph is stored in a stateless graph object which can be shared
// by multiple session instances which handle the initialization of variables related to the graph.
// Following this approach in CMSSW, a graph should be stored in a GlobalCache which can be accessed
// by sessions owned by multiple stream module copies. Instead of using only the plain graph, we
// make use of a cache struct that can be extended in the future if nedded. In addition, the graph
// is protected via std::atomic, which should not affect the performance as it is only accessed in
// the module constructor and not in the actual produce loop.
struct DeepVertexTFCache {
  DeepVertexTFCache() : graphDef(nullptr) {
  }

  std::atomic<tensorflow::GraphDef*> graphDef;
};


class DeepVertexTFJetTagsProducer : public edm::stream::EDProducer<edm::GlobalCache<DeepVertexTFCache>> {

  public:
    explicit DeepVertexTFJetTagsProducer(const edm::ParameterSet&, const DeepVertexTFCache*);
    ~DeepVertexTFJetTagsProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

    static std::unique_ptr<DeepVertexTFCache> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const DeepVertexTFCache*);

    enum InputIndexes {
      kGlobal = 0,
      kSeedingTracks = 1,
      kNeighbourTracks = 2,
    };

    enum OutputIndexes {
      kJetFlavour = 0,
      
        
    };
    
  private:
    typedef std::vector<reco::DeepFlavourTagInfo> TagInfoCollection;
    typedef reco::JetTagCollection JetTagCollection;

    void beginStream(edm::StreamID) override {}
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override {}

    const edm::EDGetTokenT< TagInfoCollection > src_;
    std::vector<std::pair<std::string,std::vector<unsigned int>>> flav_pairs_;
    std::vector<std::string> input_names_;
    std::vector<std::string> output_names_;
    std::vector<std::string> lp_names_;

    // session for TF evaluation
    tensorflow::Session* session_;
    // vector of learning phase tensors, i.e., boolean scalar tensors pointing to false
    std::vector<tensorflow::Tensor> lp_tensors_;
    // flag to evaluate model batch or jet by jet
    bool batch_eval_;
    
};


DeepVertexTFJetTagsProducer::DeepVertexTFJetTagsProducer(const edm::ParameterSet& iConfig,
  const DeepVertexTFCache* cache) :
  src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
  output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")),
  lp_names_(iConfig.getParameter<std::vector<std::string>>("lp_names")),
  session_(nullptr),
  batch_eval_(iConfig.getParameter<bool>("batch_eval"))
{
  // get threading config and build session options
  size_t nThreads = iConfig.getParameter<unsigned int>("nThreads");
  std::string singleThreadPool = iConfig.getParameter<std::string>("singleThreadPool");
  tensorflow::SessionOptions sessionOptions;
  tensorflow::setThreading(sessionOptions, nThreads, singleThreadPool);

  // create the session using the meta graph from the cache
  session_ = tensorflow::createSession(cache->graphDef, sessionOptions);

  // get output names from flav_table
  const auto & flav_pset = iConfig.getParameter<edm::ParameterSet>("flav_table");
  for (const auto flav_pair : flav_pset.tbl()) {
    const auto & flav_name = flav_pair.first;
    flav_pairs_.emplace_back(flav_name,
                             flav_pset.getParameter<std::vector<unsigned int>>(flav_name));
  }

  for (const auto & flav_pair : flav_pairs_) {
    produces<JetTagCollection>(flav_pair.first);
  }

  // flag inputs (required because of batch norm)
  // names for the learing phase placeholders (to init and set as false)
  for (size_t i = 0; i < lp_names_.size(); i++) {
    // create a bool tensor, set its value to false and store it
    tensorflow::Tensor t(tensorflow::DT_BOOL, {});
    t.scalar<bool>()() = false;
    lp_tensors_.push_back(t);
  }
}


DeepVertexTFJetTagsProducer::~DeepVertexTFJetTagsProducer()
{
  // close and delete the session
  if (session_ != nullptr) {
    tensorflow::closeSession(session_);
  }
}


void DeepVertexTFJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{

  // pfDeepVertexJetTags
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfDeepFlavourTagInfos"));
  desc.add<std::vector<std::string>>("input_names", 
    { "input_1", "input_2", "input_3",  "input_4","input_5","input_6","input_7","input_8","input_9","input_10","input_11","input_12" });
  desc.add<edm::FileInPath>("graph_path",
    edm::FileInPath("RecoBTag/TensorFlow/data/Converted_retraining.pb"));

  desc.add<std::vector<std::string>>("lp_names",
    { "globals_input_batchnorm/keras_learning_phase" });
  desc.add<std::vector<std::string>>("output_names",
    { "output_node0" }); ///SIGMOID
  {
    edm::ParameterSetDescription psd0;
    psd0.add<std::vector<unsigned int>>("probb", {0});
    psd0.add<std::vector<unsigned int>>("probc", {1});
    psd0.add<std::vector<unsigned int>>("probu", {2});
    psd0.add<std::vector<unsigned int>>("probg", {3});
    desc.add<edm::ParameterSetDescription>("flav_table", psd0);
  }

  desc.add<bool>("batch_eval", false);

  desc.add<unsigned int>("nThreads", 1);
  desc.add<std::string>("singleThreadPool", "no_threads");

  descriptions.add("pfDeepVertexJetTags", desc);
}

std::unique_ptr<DeepVertexTFCache> DeepVertexTFJetTagsProducer::initializeGlobalCache(
  const edm::ParameterSet& iConfig)
{
  // set the tensorflow log level to error
  tensorflow::setLogging("3");

  // get the pb file
  std::string pbFile = iConfig.getParameter<edm::FileInPath>("graph_path").fullPath();

  // load the graph def and save it in the cache
  DeepVertexTFCache* cache = new DeepVertexTFCache();
  cache->graphDef = tensorflow::loadGraphDef(pbFile);

  return std::unique_ptr<DeepVertexTFCache>(cache);
}

void DeepVertexTFJetTagsProducer::globalEndJob(const DeepVertexTFCache* cache)
{
  if (cache->graphDef != nullptr) {
    delete cache->graphDef;
  }
}


void DeepVertexTFJetTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);
  
  // initialize output collection
  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  for (std::size_t i=0; i < flav_pairs_.size(); i++) {
    if (!tag_infos->empty()) {
      auto jet_ref = tag_infos->begin()->jet();
      output_tags.emplace_back(std::make_unique<JetTagCollection>(
            edm::makeRefToBaseProdFrom(jet_ref, iEvent)));
    } else {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
}

////changes here:::///



  const int64_t n_jets = tag_infos->size();
  // either all jets or one per batch for the time being
  const int64_t n_batch_jets = batch_eval_ ? n_jets : 1;
  
  std::vector<tensorflow::TensorShape> input_sizes {
    {n_batch_jets, 4},         // input_1 - global jet features
 
     {n_batch_jets, 10, 21},  
    {n_batch_jets,20, 36},      // input_3 - neighbours  
    {n_batch_jets,20, 36},      // input_4 - neighbours  
    {n_batch_jets,20, 36},      // input_5 - neighbours  
    {n_batch_jets,20, 36},      // input_6 - neighbours  
    {n_batch_jets,20, 36},      // input_7 - neighbours  
    {n_batch_jets,20, 36},      // input_8 - neighbours  
    {n_batch_jets,20, 36},      // input_9 - neighbours  
    {n_batch_jets,20, 36},      // input_10 - neighbours  
    {n_batch_jets,20, 36},      // input_11 - neighbours  
    {n_batch_jets,20, 36},      // input_12 - neighbours 
/* 

    {n_batch_jets, 21, 10},     // input_2 - tracks
    {n_batch_jets, 36, 20},      // input_3 - neighbours  
    {n_batch_jets, 36, 20},      // input_4 - neighbours  
    {n_batch_jets, 36, 20},      // input_5 - neighbours  
    {n_batch_jets, 36, 20},      // input_6 - neighbours  
    {n_batch_jets, 36, 20},      // input_7 - neighbours  
    {n_batch_jets, 36, 20},      // input_8 - neighbours  
    {n_batch_jets, 36, 20},      // input_9 - neighbours  
    {n_batch_jets, 36, 20},      // input_10 - neighbours  
    {n_batch_jets, 36, 20},      // input_11 - neighbours  
    {n_batch_jets, 36, 20},      // input_12 - neighbours  */

    
  };
  
  
  // create a list of named tensors, i.e. a vector of (string, Tensor) pairs, with proper size to
  // prevent element copying that would occur via push_back's
  // the default Tensor constructor creates a scalar so this should be fine w.r.t. to memory
  tensorflow::NamedTensorList input_tensors;
  input_tensors.resize(input_sizes.size() ); //+ lp_tensors_.size());
  
  // add actual input tensors that hold physics information
  for (std::size_t i=0; i < input_sizes.size(); i++) {
    input_tensors[i] = tensorflow::NamedTensor(
      input_names_[i], tensorflow::Tensor(tensorflow::DT_FLOAT, input_sizes.at(i)));
  }
  
  // add learning-phase tensors behind them
 // for (std::size_t i=0; i < lp_tensors_.size(); i++) {
 //   input_tensors[input_sizes.size() + i] = tensorflow::NamedTensor(lp_names_[i], lp_tensors_[i]);
 // }
  
  
  
  std::size_t n_batches = n_jets/n_batch_jets; // either 1 or n_jets
  for (std::size_t batch_n=0; batch_n < n_batches; batch_n++) {

    // tensors have to be zeroed before filling per batch
    for (std::size_t i=0; i < input_sizes.size(); i++) {
      input_tensors[i].second.flat<float>().setZero();
    }

    // fill values of the input tensors
    for (std::size_t jet_bn=0; jet_bn < (std::size_t) n_batch_jets; jet_bn++) {

      // global jet index (jet_bn is the jet batch index)
      std::size_t jet_n = batch_n*n_batch_jets + jet_bn;

      // jet and other global features
      const auto & features = tag_infos->at(jet_n).features();
      
//       jvars_mean_BIS=numpy.array([65.0491047, -0.00124901765, -0.0179816545, 10.6507911, 7.12620145, 0.415715371, -0.33098271, 1.75991073, 0.351868548])
//       jvars_std_BIS=numpy.array([20.19552972, 1.20940063, 1.81650516, 3.40517432, 7.31917375, 0.6973823, 0.81846396, 21.01820818, 1.88900173])
      
      
      //add phi and mass!!
      
      input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 0) = (features.jet_features.pt);//-65.0491047)/20.19552972;
      input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 1) = (features.jet_features.eta);//+0.00124901765)/1.20940063;
      input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 2) = (features.jet_features.phi);//+0.0179816545)/1.81650516;
      input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 3) = (features.jet_features.mass);//-10.6507911)/3.40517432;
      
      
      std::cout<<input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 0)<<" "<<input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 1)<<
      " "<<input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 2)<<" "<<input_tensors.at(kGlobal).second.matrix<float>()(jet_bn, 3)<<" "<<std::endl;

      // seed features
      auto max_seed_n = std::min(features.seed_features.size(),
        (std::size_t) input_sizes.at(kSeedingTracks).dim_size(1));
      
      std::cout<<"what my size"<<" "<<max_seed_n<<" "<<std::min(features.seed_features.size(),
        (std::size_t) input_sizes.at(kSeedingTracks).dim_size(1))<<" "<<features.c_pf_features.size()<<" "<<features.seed_features.size()<<" "<<(std::size_t) input_sizes.at(kSeedingTracks).dim_size(1)<<std::endl;
      
      for (std::size_t seed_n=0; seed_n < max_seed_n; seed_n++) {
        const auto & seed_features = features.seed_features.at(seed_n);
        
        seed_tensor_filler(input_tensors.at(kSeedingTracks).second,
                           jet_bn, seed_n, seed_features);
        
        auto max_neighbour_n = std::min(features.seed_features[seed_n].seed_nearTracks.size(),
        (std::size_t) input_sizes.at(kNeighbourTracks).dim_size(1));
        
        std::cout<<"what neighbours"<<" "<<(std::size_t) input_sizes.at(kNeighbourTracks).dim_size(1)<<" "<<features.seed_features[seed_n].seed_nearTracks.size()<<std::endl;
        
 //       for (std::size_t neighbour_n=0; neighbour_n < max_neighbour_n; seed_n++) {
        neighbourTracks_tensor_filler(input_tensors.at(kNeighbourTracks+seed_n).second,
                           jet_bn, seed_n,  seed_features);
            
  //      }
      }
      
    
    std::cout<<"seedingtracks"<<std::endl;
    
    for (unsigned int f=0; f<10; f++){ 
        for (unsigned int fj=0; fj<21; fj++){ 
            std::cout<<input_tensors.at(kSeedingTracks).second.tensor<float, 3>()(jet_bn, f, fj)<<" ";}
        std::cout<<" "<<" "<<std::endl;  }
     
     
    std::cout<<"first neighbour"<<std::endl;
    
    for (unsigned int f=0; f<20; f++){ 
        for (unsigned int fj=0; fj<36; fj++){ 
            std::cout<<input_tensors.at(kNeighbourTracks).second.tensor<float, 3>()(jet_bn, f, fj)<<" ";}
        std::cout<<" "<<" "<<std::endl;  }


    std::cout<<"second neighbour"<<std::endl;
    
    for (unsigned int f=0; f<20; f++){ 
        for (unsigned int fj=0; fj<36; fj++){ 
            std::cout<<input_tensors.at(kNeighbourTracks+1).second.tensor<float, 3>()(jet_bn, f, fj)<<" ";}
        std::cout<<" "<<" "<<std::endl;  }
        
    }////different

    

    

    // run the session
    std::vector<tensorflow::Tensor> outputs;
//    std::cout <<"Input size" <<  input_tensors.size() << std::endl;
    tensorflow::run(session_, input_tensors, output_names_, &outputs);

    // set output values for flavour probs
    for (std::size_t jet_bn=0; jet_bn < (std::size_t) n_batch_jets; jet_bn++) {

      // global jet index (jet_bn is the jet batch index)
      std::size_t jet_n = batch_n*n_batch_jets + jet_bn;

      const auto & jet_ref = tag_infos->at(jet_n).jet();
      for (std::size_t flav_n=0; flav_n < flav_pairs_.size(); flav_n++) {
        const auto & flav_pair = flav_pairs_.at(flav_n);
        float o_sum = 0.;
        for (const unsigned int & ind : flav_pair.second) {
          o_sum += outputs.at(kJetFlavour).matrix<float>()(jet_bn, ind);
        }
        (*(output_tags.at(flav_n)))[jet_ref] = o_sum;
      }
    }
  }  
  
  

/////



  for (std::size_t i=0; i < flav_pairs_.size(); i++) {
    iEvent.put(std::move(output_tags[i]), flav_pairs_.at(i).first);
  }

}




//define this as a plug-in
DEFINE_FWK_MODULE(DeepVertexTFJetTagsProducer);


