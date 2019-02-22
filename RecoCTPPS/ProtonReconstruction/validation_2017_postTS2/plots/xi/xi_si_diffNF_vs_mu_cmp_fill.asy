import root;
import pad_layout;
include "../defaults.asy";

string topDir = "../../data_eos/";

string datasets[] = def_datasets;

string stream = def_stream;

string alignment = def_alignment;

string cols[], c_labels[], c_si_rps[];
cols.push("arm0"); c_labels.push("sector 45 (L)"); c_si_rps.push("rp23");
cols.push("arm1"); c_labels.push("sector 56 (R)"); c_si_rps.push("rp123");

//xTicksDef = LeftTicks(0.05, 0.01);
//yTicksDef = RightTicks(0.05, 0.01);

//----------------------------------------------------------------------------------------------------

NewPad(false);
label("\vbox{\hbox{stream: " + stream + "}\hbox{alignment: " + replace(alignment, "_", "\_") + "}}");

for (int ci : cols.keys)
	NewPadLabel(c_labels[ci]);

for (int dsi : datasets.keys)
{
	NewRow();

	NewPadLabel(replace(datasets[dsi], "_", "\_"));

	for (int ci : cols.keys)
	{
		NewPad("$\xi_{\rm multi}$", "$\xi_{\rm single,F} - \xi_{\rm single,N}$");

		string f = topDir + datasets[dsi] + "/" + stream + "/alignment_" + alignment + "/output.root";
		string on = "armCorrelationPlots/" + cols[ci] + "/p_xi_si_diffNF_vs_xi_mu";
		
		RootObject hist = RootGetObject(f, on, error=false);

		if (!hist.valid)
			continue;
		
		draw(hist, "eb,d0", red);

		limits((0., -0.02), (0.25, +0.02), Crop);
	}
}

GShipout("xi_si_diffNF_vs_mu_cmp_fill", hSkip=1mm, vSkip=0mm);
