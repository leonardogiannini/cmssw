/****************************************************************************
 *
 * This is a part of CTPPS offline software.
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/

#ifndef DataFormats_ProtonReco_ProtonTrack_h
#define DataFormats_ProtonReco_ProtonTrack_h

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include <set>

/**
 * FIXME make use of the general reco::Candidate object, with appropriate set'ters and get'ters
 */

namespace reco
{
  class ProtonTrack
  {
    public:
      ProtonTrack() :
        xi_( 0. ), xi_unc_( 0. ),
        isValid_( false ) {}
      ProtonTrack( const Local3DPoint& vtx, const Local3DVector& dir, float xi, float xi_unc=0. ) :
        vertex_( vtx ), direction_( dir ), xi_( xi ), xi_unc_( xi_unc ),
        isValid_( true ) {}
      ~ProtonTrack() {}

      // in mm
      void setVertex( const Local3DPoint& vtx ) { vertex_ = vtx; }
      const Local3DPoint& vertex() const { return vertex_; }

      // momentum in GeV
      void setDirection( const Local3DVector& dir ) { direction_ = dir; }
      const Local3DVector& direction() const { return direction_; }

      void setXi( float xi ) { xi_ = xi; }
      float xi() const { return xi_; }
      float xiErr() const { return xi_unc_; }

      void setValid( bool valid=true ) { isValid_ = valid; }
      bool valid() const { return isValid_; }

      // TODO: add proper getters, setters
      enum { rmSingleRP, rmMultiRP } method;
      enum { sector45, sector56 } lhcSector;

      bool singleRP() const { return method==rmSingleRP; }
      bool multiRP() const { return method==rmMultiRP; }
      bool isSector45() const { return lhcSector==sector45; }
      bool isSector56() const { return lhcSector==sector56; } 
      
        // static const int kDetMask            = 0xF;
        // static const int kSubdetMask         = 0x7;
        // static const int kDetOffset          = 28;
        // static const int kSubdetOffset       = 25;
        // static const int startArmBit         = 24;
        // static const int maskArm             = 0x1;
        // static const int startStationBit     = 22;
        // static const int maskStation         = 0x3;
        // static const int startRPBit          = 19;
        // static const int maskRP              = 0x7;
        
      bool pixelFlag() const
      {                 
          for(auto rp : contributingRPIds)  {if (((rp>>25)&0x7)==4) return true;}
          return false; 
      }
      
      bool armFlag() const
      {                 
          for(auto rp : contributingRPIds)  {if (((rp>>24)&0x1)==1) return true;}
          return false; 
      }
      
      bool statFlag() const
      {                 
          for(auto rp : contributingRPIds)  {if (((rp>>22)&0x3)==2) return true;}
          return false; 
      }
      
      bool rpFlag() const
      {                 
          for(auto rp : contributingRPIds)  {if (((rp>>19)&0x7)==3) return true;}
          return false; 
      }
      
      int rpArmStatId() const
      {   
          int multiplier=1;
          int res=0;
          for(auto rp : contributingRPIds) 
          {
              res=res+multiplier*(100*((rp>>24)&0x1) + 10*((rp>>22)&0x3) + ((rp>>19)&0x7));
              multiplier=multiplier*1000;
          } 
          return res;
      }
      
      
      

      unsigned int fitNDF;
      double fitChiSq;

      std::set<unsigned int> contributingRPIds;

      double halfCrossingAngleSector45=0., halfCrossingAngleSector56=0.;

    private:

      // TODO: describe, mention CMS coordinate notation
      Local3DPoint vertex_;
      Local3DVector direction_;

      // TODO: describe
      float xi_;
      float xi_unc_;

      // TODO: rename to fit valid?
      bool isValid_;
  };
}

#endif
