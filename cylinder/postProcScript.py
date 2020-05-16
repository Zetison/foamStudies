# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import os
from math import pi
import numpy as np
#### disable automatic camera reset on 'Show'
dir = os.getcwd()
# load plugin
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
cylinderfoam = GetActiveSource()
u_inf = 1
d = 1
L = 2*pi*d
denominator = 0.5*u_inf**2*L*d # (rho is canceled by the rho in the kinematik pressure)

# destroy cylinderfoam
Delete(cylinderfoam)
del cylinderfoam

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'OpenFOAMReader'
cylinderfoam = OpenFOAMReader(FileName=dir+'/cylinder.foam')
cylinderfoam.MeshRegions = ['internalMesh']
cylinderfoam.CellArrays = ['U', 'p']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1934, 1154]

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

LoadPalette(paletteName='BlackBackground')

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.7)

# set active view
SetActiveView(None)

# split cell
layout1.SplitVertical(2, 0.5)

# show data in view
cylinderfoamDisplay = Show(cylinderfoam, renderView1, 'UnstructuredGridRepresentation')

# create a new 'Gradient Of Unstructured DataSet'
gradientOfUnstructuredDataSet1 = GradientOfUnstructuredDataSet(Input=cylinderfoam)
gradientOfUnstructuredDataSet1.ScalarArray = ['POINTS', 'U']
gradientOfUnstructuredDataSet1.ComputeVorticity = 1

# show data in view
gradientOfUnstructuredDataSet1Display = Show(gradientOfUnstructuredDataSet1, renderView1, 'UnstructuredGridRepresentation')

cylinderfoamDisplay.SetScalarBarVisibility(renderView1, False)
Hide(cylinderfoam, renderView1)

# Properties modified on gradientOfUnstructuredDataSet1
gradientOfUnstructuredDataSet1.ResultArrayName = 'Vorticity'

# set scalar coloring
ColorBy(gradientOfUnstructuredDataSet1Display, ('POINTS', 'Vorticity', 'Magnitude'))

# show color bar/color legend
gradientOfUnstructuredDataSet1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('Vorticity')

# get color transfer function/color map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
vorticityLUT.ApplyPreset('SINTEF1', True)
vorticityPWF.ApplyPreset('SINTEF1', True)

# reset view to fit data
# Rescale transfer function
vorticityPWF.RescaleTransferFunction(0.0, 10.0)

# Rescale transfer function
vorticityLUT.AutomaticRescaleRangeMode = "Never"
vorticityLUT.RescaleOnVisibilityChange = 0
vorticityLUT.RescaleTransferFunction(0.0, 10.0)

# get color legend/bar for vorticityLUT in view renderView1
vorticityLUTColorBar = GetScalarBar(vorticityLUT, renderView1)

# change scalar bar placement
vorticityLUTColorBar.Orientation = 'Vertical'
vorticityLUTColorBar.WindowLocation = 'AnyLocation'
vorticityLUTColorBar.Position = [0.017238868081564718, 0.0494947987803247]
vorticityLUTColorBar.ScalarBarLength = 0.33000000000000007

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.7899, 0.0, 73.31538059794683]
renderView1.CameraFocalPoint = [5.7899, 0.0, 3.1415927410125732]
renderView1.CameraParallelScale = 8.472852946250722
renderView1.Update()
SetActiveSource(None)

if True: # plot SINTEF logo
  # create a new 'Logo'
  logo1 = Logo()
  SINTEF_white = CreateTexture("/home/zetison/OneDrive/work/graphics/logos/SINTEF_white.png")
  logo1.Texture = SINTEF_white
  logo1Display = Show(logo1, renderView1, 'LogoSourceRepresentation')
  logo1Display.Position = [0.840060240963855, 0.0013581890812250329]
  logo1Display.Interactivity = 0

#############################################################################################################
## Create 2D plots of lift and drag
# create a new 'OpenFOAMReader'
cylinderfoam_1 = OpenFOAMReader(FileName=dir+'/cylinder.foam')
cylinderfoam_1.MeshRegions = ['internalMesh']
cylinderfoam_1.CellArrays = ['U', 'p']

# Properties modified on cylinderfoam_1
cylinderfoam_1.MeshRegions = ['cylinder']

gradientOfUnstructuredDataSet2 = GradientOfUnstructuredDataSet(Input=cylinderfoam_1)
gradientOfUnstructuredDataSet2.ScalarArray = ['CELLS', 'U']
gradientOfUnstructuredDataSet2.ResultArrayName = 'dU'

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(Input=gradientOfUnstructuredDataSet2)
generateSurfaceNormals1.ComputeCellNormals = 1

# create a new 'Calculator'
calculator1 = Calculator(Input=generateSurfaceNormals1)
calculator1.AttributeType = 'Cell Data'

# Properties modified on calculator1
calculator1.ResultArrayName = 'drag'

# Properties modified on calculator1
calculator1.Function = '(p*Normals_X + 2*0.01*dU_0*Normals_X + 0.01*(dU_1+dU_3)*Normals_Y)/'+str(denominator)
#calculator1.Function = '(p*Normals_X + 2*dU_0*Normals_X + (dU_1+dU_3)*Normals_Y)/'+str(denominator)
#calculator1.Function = '(0.01*dU_0*Normals_X + 0.01/2*(dU_3+dU_1)*Normals_Y)/'+str(denominator)
#calculator1.Function = '(0.01*dU_8*Normals_X)/'+str(denominator)

# create a new 'Calculator'
calculator2 = Calculator(Input=calculator1)
calculator2.AttributeType = 'Cell Data'

# Properties modified on calculator2
calculator2.ResultArrayName = 'lift'

# Properties modified on calculator2
#calculator2.Function = 'p*Normals_Y/'+str(denominator)
calculator2.Function = '(p*Normals_Y + 2*0.01*dU_4*Normals_Y + 0.01*(dU_1+dU_3)*Normals_X)/'+str(denominator)

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=calculator2)

CreateLayout('Layout #2')

# get layout
layout2 = GetLayoutByName("Layout #2")

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')

# show data in view
integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')
spreadSheetView1.FieldAssociation = 'Cell Data'

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout2, hint=0)

SelectIDs(IDs=[-1, 0], FieldType=0, ContainingCells=0)

layout1 = GetLayoutByName("Layout #1")
# set active source
SetActiveSource(integrateVariables1)

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=quartileChartView1, layout=layout1, hint=1)

# uncomment following to set a specific view size
# quartileChartView1.ViewSize = [400, 400]
plotSelectionOverTime1 = PlotSelectionOverTime(Input=integrateVariables1, Selection=None)
# show data in view
plotSelectionOverTime1Display = Show(plotSelectionOverTime1, quartileChartView1, 'QuartileChartRepresentation')

# update the view to ensure updated data information
quartileChartView1.Update()

# Properties modified on plotSelectionOverTime1Display
plotSelectionOverTime1Display.SeriesVisibility = ['drag (stats)']

# Properties modified on quartileChartView1
quartileChartView1.LeftAxisUseCustomRange = 1

# Properties modified on quartileChartView1
quartileChartView1.LeftAxisRangeMaximum = 1.4

# Properties modified on quartileChartView1
quartileChartView1.LeftAxisRangeMinimum = 1.0

# Properties modified on quartileChartView1
quartileChartView1.ShowLegend = 0

# Properties modified on quartileChartView1
quartileChartView1.LeftAxisTitle = 'Drag coefficient'

# Properties modified on quartileChartView1
quartileChartView1.BottomAxisTitle = 'Time [s]'

# Properties modified on plotSelectionOverTime2Display
plotSelectionOverTime1Display.SeriesColor = ['drag (stats)', '0', '0', '0']

# create a new 'CSV Reader'
coefficientcsv = CSVReader(FileName=['/home/zetison/OpenFOAM/OpenFOAM-v1912/run/cylinder/postProcessing/forces/0/coefficient.csv'])

# show data in view
coefficientcsvDisplay_1 = Show(coefficientcsv, quartileChartView1, 'QuartileChartRepresentation')
coefficientcsvDisplay_1.SeriesVisibility = ['Cd']
coefficientcsvDisplay_1.SeriesColor = ['Cd', '1', '0', '0']

# Create a new 'Quartile Chart View'
quartileChartView2 = CreateView('QuartileChartView')
# uncomment following to set a specific view size
# quartileChartView2.ViewSize = [400, 400]

plotSelectionOverTime2 = PlotSelectionOverTime(Input=integrateVariables1, Selection=None)
# show data in view
plotSelectionOverTime2Display = Show(plotSelectionOverTime2, quartileChartView2, 'QuartileChartRepresentation')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=quartileChartView2, layout=layout1, hint=2)

# update the view to ensure updated data information
quartileChartView1.Update()

# update the view to ensure updated data information
quartileChartView2.Update()

# Properties modified on plotSelectionOverTime2Display
plotSelectionOverTime2Display.SeriesVisibility = ['lift (stats)']

# Properties modified on quartileChartView1
quartileChartView2.LeftAxisUseCustomRange = 1

# Properties modified on quartileChartView1
quartileChartView2.LeftAxisRangeMaximum = 0.4

# Properties modified on quartileChartView1
quartileChartView2.LeftAxisRangeMinimum = -0.4

# Properties modified on quartileChartView2
quartileChartView2.ShowLegend = 0

# Properties modified on quartileChartView2
quartileChartView2.BottomAxisTitle = 'Time [s]'

# Properties modified on quartileChartView2
quartileChartView2.LeftAxisTitle = 'Lift coefficient'

# Properties modified on plotSelectionOverTime2Display
plotSelectionOverTime2Display.SeriesColor = ['lift (stats)', '0', '0', '0']

coefficientcsvDisplay_2 = Show(coefficientcsv, quartileChartView2, 'QuartileChartRepresentation')
coefficientcsvDisplay_2.SeriesVisibility = ['Cl']
coefficientcsvDisplay_2.SeriesColor = ['Cl', '1', '0', '0']

# set active view
SetActiveView(spreadSheetView1)
SelectIDs(IDs=[-1, 0], FieldType=0, ContainingCells=0)
ExportView('/home/zetison/OpenFOAM/results/cylinder.csv', view=spreadSheetView1)
####################################################################################################
# resize frame
SetActiveView(renderView1)

RenderAllViews()
if False:
	# get animation scene
	SaveAnimation('/home/zetison/Videos/animation.ogv', layout1, 
	    FontScaling='Scale fonts proportionally',
	    OverrideColorPalette='',
	    StereoMode='No change',
	    TransparentBackground=0,
      SaveAllViews=1,
	    ImageQuality=100,
	    FrameRate=15,
      ImageResolution=[1920, 1080],
      SeparatorWidth=0,
      SeparatorColor=[1.0, 1.0, 1.0],
			Quality=2) #,
#      FrameWindow=[0,10],

#### uncomment the following to render all views
# alternatively, if you want to write images, you can use SaveScreenshot(...).
