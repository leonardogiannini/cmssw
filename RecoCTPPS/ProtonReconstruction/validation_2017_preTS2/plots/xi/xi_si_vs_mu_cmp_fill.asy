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

xTicksDef = LeftTicks(0.05, 0.01);
yTicksDef = RightTicks(0.05, 0.01);

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
		NewPad("$\xi_{\rm single}$", "$\xi_{\rm multi}$");

		string f = topDir + datasets[dsi] + "/" + stream + "/alignment_" + alignment + "/output.root";
		string on = "singleMultiCorrelationPlots/si_" + c_si_rps[ci]  + "_mu_" + cols[ci] + "/h2_xi_mu_vs_xi_si";
		
		RootObject hist = RootGetObject(f, on, error=false);

		if (!hist.valid)
			continue;
		
		draw(hist);

		draw((0, 0)--(0.2, 0.2), dashed);

		limits((0., 0.), (0.2, 0.2), Crop);
	}
}

GShipout("xi_si_vs_mu_cmp_fill", hSkip=1mm, vSkip=0mm);
