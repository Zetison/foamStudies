#### import the simple module from the paraview
from paraview.simple import *
import os
import numpy as np
import getpass
#### disable automatic camera reset on 'Show'
dir = os.getcwd()
# load plugin
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

viewSize = [1920, 1080]
axisLabelFontSize=13
fontSize=18
staticCylinder=True
username = getpass.getuser()
home = '/home/'+username+'/'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'OpenFOAMReader'
SED_MODEL_NAMEfoam = OpenFOAMReader(FileName=dir+'/SED_MODEL_NAME.foam')
SED_MODEL_NAMEfoam.CaseType = 'Decomposed Case'
SED_MODEL_NAMEfoam.MeshRegions = ['internalMesh']
SED_MODEL_NAMEfoam.CellArrays = ['vorticity']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

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
SED_MODEL_NAMEfoamDisplay = Show(SED_MODEL_NAMEfoam, renderView1, 'UnstructuredGridRepresentation')

# set scalar coloring
ColorBy(SED_MODEL_NAMEfoamDisplay, ('POINTS', 'vorticity', 'Magnitude'))

# show color bar/color legend
SED_MODEL_NAMEfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('vorticity')

# get color transfer function/color map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('vorticity')

with open(home+"kode/colormaps/SINTEF1.xml", "r") as f:
      data = f.read()
      vorticityLUT.ApplyColorMap(data)
#vorticityLUT.ApplyPreset('SINTEF1', True)
#vorticityPWF.ApplyPreset('SINTEF1', True)

# reset view to fit data
# Rescale transfer function
vorticityPWF.RescaleTransferFunction(0.0, 10.0)

# Rescale transfer function
vorticityLUT.AutomaticRescaleRangeMode = "Never"
vorticityLUT.RescaleOnVisibilityChange = 0
if staticCylinder:
	vorticityLUT.RescaleTransferFunction(0.0, 10.0)
else:
	vorticityLUT.RescaleTransferFunction(0.0, 20.0)

# get color legend/bar for vorticityLUT in view renderView1
vorticityLUTColorBar = GetScalarBar(vorticityLUT, renderView1)

# change scalar bar placement
vorticityLUTColorBar.Orientation = 'Vertical'
vorticityLUTColorBar.WindowLocation = 'LowerLeftCorner'
vorticityLUTColorBar.Title = 'Vorticity'
vorticityLUTColorBar.ComponentTitle = 'magnitude'
#vorticityLUTColorBar.TitleFontSize = fontSize
#vorticityLUTColorBar.LabelFontSize = fontSize
#vorticityLUTColorBar.ScalarBarThickness = 10
#vorticityLUTColorBar.ScalarBarLength = 0.3

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.7899, 0.0, 73.31]
renderView1.CameraFocalPoint = [5.7899, 0.0, 3.1415]
renderView1.CameraParallelScale = 8.4728

if True: # plot SINTEF logo
	# create a new 'Logo'
	logo1 = Logo()
	logo1.Texture = CreateTexture(home+'OneDrive/work/graphics/logos/SINTEF_white.png')
	logo1Display = Show(logo1, renderView1, 'LogoSourceRepresentation')
	logo1Display.Position = [0.86, 0.0]
	logo1Display.Interactivity = 0

#############################################################################################################
## Add volleyball image
if not(staticCylinder):
	plane1 = Plane()
	D = 1.0
	plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')
	plane1.Origin = [-D/2, -D/2, 10.0]
	plane1.Point1 = [D/2, -D/2, 10.0]
	plane1.Point2 = [-D/2, D/2, 10.0]
	transform1 = Transform(Input=plane1)
	volleyball = CreateTexture(home+"OneDrive/work/paraview/sources/volleyball.png")
	transform1Display = Show(transform1, renderView1, 'GeometryRepresentation')
	Hide(plane1,renderView1)
	transform1Display.Texture = volleyball 
	transform1TransformRotationTrack = GetAnimationTrack('Rotation', index=2, proxy=transform1.Transform)
#	keyFrame10769 = CompositeKeyFrame()
#	keyFrame10770.KeyTime = 1.0
#	transform1TransformRotationTrack.KeyFrames = [keyFrame10769, keyFrame10770]
#############################################################################################################
## Create 2D plots of lift and drag

# create a new 'CSV Reader'
coefficientcsv = CSVReader(FileName=[dir+'/postProcessing/forceCoeffs/0/coefficient.csv'])

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=quartileChartView1, layout=layout1, hint=1)

# Properties modified on quartileChartView1
if staticCylinder:
	quartileChartView1.LeftAxisUseCustomRange = 1
	quartileChartView1.LeftAxisRangeMaximum = 1.4
	quartileChartView1.LeftAxisRangeMinimum = 1.0
else:
	quartileChartView1.LeftAxisUseCustomRange = 1
	quartileChartView1.LeftAxisRangeMaximum = 3.0
	quartileChartView1.LeftAxisRangeMinimum = 0.0

quartileChartView1.ShowLegend = 0
quartileChartView1.LeftAxisTitle = 'Drag coefficient'
quartileChartView1.BottomAxisTitle = 'Time [s]'
#quartileChartView1.LeftAxisTitleFontSize = fontSize
#quartileChartView1.BottomAxisTitleFontSize = fontSize
quartileChartView1.LeftAxisLabelFontSize = axisLabelFontSize
quartileChartView1.BottomAxisLabelFontSize = axisLabelFontSize

# show data in view
coefficientcsvDisplay_1 = Show(coefficientcsv, quartileChartView1, 'QuartileChartRepresentation')
coefficientcsvDisplay_1.SeriesVisibility = ['Cd']
coefficientcsvDisplay_1.SeriesColor = ['Cd', '0', '0', '0']
quartileChartView1.Update()

# Create a new 'Quartile Chart View'
quartileChartView2 = CreateView('QuartileChartView')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=quartileChartView2, layout=layout1, hint=2)

# Properties modified on quartileChartView1
if staticCylinder:
	quartileChartView2.LeftAxisUseCustomRange = 1
	quartileChartView2.LeftAxisRangeMaximum = 0.4
	quartileChartView2.LeftAxisRangeMinimum = -0.4

quartileChartView2.ShowLegend = 0
quartileChartView2.BottomAxisTitle = 'Time [s]'
quartileChartView2.LeftAxisTitle = 'Lift coefficient'
#quartileChartView2.LeftAxisTitleFontSize = fontSize
#quartileChartView2.BottomAxisTitleFontSize = fontSize
quartileChartView2.LeftAxisLabelFontSize = axisLabelFontSize
quartileChartView2.BottomAxisLabelFontSize = axisLabelFontSize

# show data in view
coefficientcsvDisplay_2 = Show(coefficientcsv, quartileChartView2, 'QuartileChartRepresentation')
coefficientcsvDisplay_2.SeriesVisibility = ['Cl']
coefficientcsvDisplay_2.SeriesColor = ['Cl', '0', '0', '0']
quartileChartView2.Update()
SetActiveView(quartileChartView2)

####################################################################################################
# Export visualization
renderView1.Update()
SetActiveSource(None)
Render()

if True:
	## Save snapshots
	for t in [100, 200, 1000]:
		timeKeeper1 = GetTimeKeeper()
		animationScene1 = GetAnimationScene()
		animationScene1.AnimationTime = t
		timeKeeper1.Time = t
		renderView1.ViewSize = viewSize
		SaveScreenshot(dir+'/screenshot_t'+str(t)+'.png', layout1, 
				FontScaling='Scale fonts proportionally',
				OverrideColorPalette='',
				StereoMode='No change',
				TransparentBackground=1,
				ImageResolution=viewSize,
				ImageQuality=100)

if False:
	# Generate movie
	SaveAnimation(dir+'/animation.ogv', layout1, 
			FontScaling='Scale fonts proportionally',
			OverrideColorPalette='',
			StereoMode='No change',
			TransparentBackground=0, 
			SaveAllViews=1,
			ImageQuality=100,
			FrameRate=25,
			ImageResolution=viewSize,
			SeparatorWidth=0,
			SeparatorColor=[1.0, 1.0, 1.0],
			Quality=2) #,
#      FrameWindow=[0,10],
