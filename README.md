# MouseAtlas
Tool to interactively create a Mouse Brain Atlas. 

You need Matlab (tested on v. 2016) and SPM12.

## Installation
Copy the files to a location of your choice.
Make sure that you included your SPM12 path in Matlab.
Open the application by 'ms_gui_mouseAtlas'

## Creating an atlas
- On the left, there is a tree selection (starting with 'root'). From there you can explore the atlas.
- To 'select' a region right click its name in the tree selection. Selected regions are getting a bold green font.
- For a visual presentation you can click 'Update'. The axis on the right is updated. The color of a region of interest (ROI) is the same as indicated in the selection tree.
- Click 'Create Atlas'. After creation you can inspect the result (a figure opens automatically). In total, the app will create 3 files:
  * the annotation atlas (matching the Allen avg template)
  * a txt file linking the numbers with the names of regions
  * the Allen atlas average template
  * an annotation file in paxinos space (suffix '_inPax')
