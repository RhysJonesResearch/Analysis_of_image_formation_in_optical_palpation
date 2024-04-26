Submitted for consideration as: Jones, R., Fang, Q. & Kennedy, B. F. (2024). Analysis of image formation in optical palpation. Journal of Biophotonics.

Summary: Code and finite element models to model optical palpation, an emerging elastography technique to image mechanical stress in 2D across a surface, with applications in intraoperative cancer detection. The analysis investigates trends in detectability, feature resolution, and contrast ratio and is validated experiments on silicone phantoms.

To generate FEA data, the CompressionOpticalPalpation.cae file was used in Abaqus, with the model 'Mesh17-Final'. The code to submit all combinations of variables as jobs is available as a macro in the file 'abaqusMacros.py' with function title 'Final_Mesh17_DeleteFilesAfterJob'.

The analysis was undertaken using the jupyter notebook files begining with 'Final_Import_and_Analysis-UseBoundaryConditions'.

'S_Dual_OP_from_ThorIm' is a matlab script that takes one through importing OCT images for Thorlabs 'ThorIm' imaging software to matlab and then implementing optical palpation in either a dual arm or common path implementation depending on the layer thickness used. 'Init_Config_231104_Inc_5_0p5mm3' is a config file that needs to be filled out and run before the 'S_Dual_OP_from_ThorIm' script. It is recommended to use the 'Plot_BScans' script to visualise the OCT data before beginning other processing, it also contains the ocde to remove bright lines in BScans that are constant throughout an entire image.

The experimental validation data and analysis is presented in 'Experimental_Validation_Plots'.

For more detailed information about the use of computational optical palpation visit https://github.com/philipwijesinghe/computational-optical-palpation.

Requirements
- matlab (tested with 2016a), including 3rd party toolboxes:
   - gridfit
   - Inpaint_nans
- Jupyter Notebook, python (tested with 3.6)
- Abaqus (tested with 6.13)
