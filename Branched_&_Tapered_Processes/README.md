FILE DESCRIPTION:

This folder includes scripts that are run stand-alone (i.e., simulations not controlled by a higher-level script).
Comments describe individual parameters and the computations carried out.

The following scripts are used to set parameters for individual simulation runs:<br>
(Parameters are set by hand in the script by the user.)<br>
* setParametersTaperProc.m <br>
	Parameter generation for a tapered process with a 2-D diffusion model.<br>
* setParametersTaperProc3D.m<br>
 	Parameter generation for a tapered process with a 3-D diffusion model<br>
	(for use in conjunction with ryanodine receptor placement at various circumferential angles).<br>
* setParametersBranchProc.m<br>
	Parameter generation for a process including a branch.<br>

The following scripts run simulations:<br>
* ProcCaDynamicsTaper.m<br>
	Simulates a tapered process.<br>
* ProcCaDynamicsTaper_Ectrl.m<br>
	Simulates a tapered process, in which the resting ER calcium level is controlled with a virtual Ca source.<br>
* ProcCaDynamicsTaper3D_ranRyR.m<br>
	Simulates a tapered process with 3-D diffusion model and ryanodine receptors placed randomly around circumference.<br>
* ProcCaDynamicsMultiSeg.m<br>
	Simulates a process with multiple segments; used to evaluate branching processes.<br>
