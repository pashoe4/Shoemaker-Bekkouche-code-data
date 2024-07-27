FILE DESCRIPTION:

This folder includes scripts that are run stand-alone (i.e., simulations not controlled by a higher-level script).
Comments describe individual parameters and the computations carried out.

The following scripts are used to set parameters for individual simulation runs.
Parameters are set by hand in the script by the user.

setParametersTaperProc.m
	Parameter generation for a tapered process with a 2-D diffusion model.
setParametersTaperProc3D.m
 	Parameter generation for a tapered process with a 3-D diffusion model
	(for use in conjunction with ryanodine receptor placement at various circumferential angles).
setParametersBranchProc.m
	Parameter generation for a process including a branch.

The following scripts run simulations.

ProcCaDynamicsTaper.m
	Simulates a tapered process.
ProcCaDynamicsTaper_Ectrl.m
	Simulates a tapered process, in which the resting ER calcium level is controlled with a virtual Ca source.
ProcCaDynamicsTaper3D_ranRyR.m
	Simulates a tapered process with 3-D diffusion model and ryanodine receptors placed randomly around circumference.
ProcCaDynamicsMultiSeg.m
	Simulates a process with multiple segments; used to evaluate branching processes.