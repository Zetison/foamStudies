#### import the simple module from the paraview
from paraview.simple import *
import os
from math import pi
import numpy as np
import getpass
#### disable automatic camera reset on 'Show'
dir = os.getcwd()
# load plugin
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

username = getpass.getuser()
home = '/home/'+username+'/'
# get active source.
cylinderfoam = GetActiveSource()

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
cylinderfoam.CellArrays = ['vorticity']

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

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

# set scalar coloring
ColorBy(cylinderfoamDisplay, ('POINTS', 'vorticity', 'Magnitude'))

# show color bar/color legend
cylinderfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('vorticity')

# get color transfer function/color map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('vorticity')

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
vorticityLUTColorBar.WindowLocation = 'LowerLeftCorner'
vorticityLUTColorBar.Title = 'Velocity'
vorticityLUTColorBar.ComponentTitle = 'magnitude'

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.7899, 0.0, 73.31]
renderView1.CameraFocalPoint = [5.7899, 0.0, 3.1415]
renderView1.CameraParallelScale = 8.4728
renderView1.Update()
SetActiveSource(None)

if True: # plot SINTEF logo
  # create a new 'Logo'
  logo1 = Logo()
  SINTEF_white = CreateTexture(home+'OneDrive/work/graphics/logos/SINTEF_white.png')
  logo1.Texture = SINTEF_white
  logo1Display = Show(logo1, renderView1, 'LogoSourceRepresentation')
  logo1Display.Position = [0.84, 0.0]
  logo1Display.Interactivity = 0

#############################################################################################################
## Add volleyball image
if True:
	plane1 = Plane()
	plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')
	plane1.Origin = [-0.5, -0.5, 10.0]
	plane1.Point1 = [0.5, -0.5, 10.0]
	plane1.Point2 = [-0.5, 0.5, 10.0]
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
staticCylinder=True
if staticCylinder:
	quartileChartView1.LeftAxisUseCustomRange = 1
	quartileChartView1.LeftAxisRangeMaximum = 1.4
	quartileChartView1.LeftAxisRangeMinimum = 1.0
quartileChartView1.ShowLegend = 0
quartileChartView1.LeftAxisTitle = 'Drag coefficient'
quartileChartView1.BottomAxisTitle = 'Time [s]'

# show data in view
coefficientcsvDisplay_1 = Show(coefficientcsv, quartileChartView1, 'QuartileChartRepresentation')
coefficientcsvDisplay_1.SeriesVisibility = ['Cd']
coefficientcsvDisplay_1.SeriesColor = ['Cd', '0', '0', '0']

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

# show data in view
coefficientcsvDisplay_2 = Show(coefficientcsv, quartileChartView2, 'QuartileChartRepresentation')
coefficientcsvDisplay_2.SeriesVisibility = ['Cl']
coefficientcsvDisplay_2.SeriesColor = ['Cl', '0', '0', '0']

####################################################################################################
# Export visualization
SetActiveView(renderView1)

RenderAllViews()
if False:
	# get animation scene
	SaveAnimation(dir+'/animation.ogv', layout1, 
	    FontScaling='Scale fonts proportionally',
	    OverrideColorPalette='',
	    StereoMode='No change',
	    TransparentBackground=0,
      SaveAllViews=1,
	    ImageQuality=100,
	    FrameRate=25,
      ImageResolution=[1920, 1080],
      SeparatorWidth=0,
      SeparatorColor=[1.0, 1.0, 1.0],
			Quality=2) #,
#      FrameWindow=[0,10],

