//Uses the generated MoleculeArchive peak tracker information and displays ROI circles for all tracked traces per slice.
//This enables the user to check the tracking results.
// Use Image>Overlay>Labels to show the UID numbering on the Image

#@ ImagePlus image
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import java.awt.Color;
import ij.gui.PointRoi;
import ij.gui.Overlay;
import java.util.Random;

Overlay overlay = new Overlay();
Random ran = new Random();

archive.getMoleculeUIDs().stream().forEach({ UID ->
   Color color = new Color(ran.nextFloat(), ran.nextFloat(), ran.nextFloat())
   table = archive.get(UID).getDataTable()
   for (int i=0; i< table.getRowCount() ; i++ ) {
      PointRoi peakRoi = new PointRoi(table.getValue("x", i) + 0.5, table.getValue("y", i) + 0.5)
      //Have to select Use Names as Labels in Image>Overlay>Labels...
      peakRoi.setName(UID.substring(0, 5))
      peakRoi.setStrokeColor(color)
      //0-4 as options...
      peakRoi.setSize(4)
      //(0=hybrid, 1=crosshair, 2=dot, 3=circle)
      peakRoi.setPointType(2)
      peakRoi.setPosition((int)table.getValue("slice", i))
      overlay.add(peakRoi)
      }
});
image.setOverlay(overlay);

//To show traces specifically per tag use the following script
#@ ImagePlus image
#@ MoleculeArchive archive
#@ String(label="Tags (comma separated list)") Tags

import de.mpg.biochem.mars.molecule.*;
import de.mpg.biochem.mars.table.*;
import java.awt.Color;
import ij.gui.PointRoi;
import ij.gui.Overlay;
import java.util.Random;

String[] tagList = Tags.split(",")

for (int i=0; i<tagList.length; i++) {
   tagList[i] = tagList[i].trim();
}

ArrayList<String> moleculesWithTags = (ArrayList<String>)archive.getMoleculeUIDs().stream().filter{UID ->
                    int tagCount = 0
                    for (int i=0; i<tagList.length; i++) {
                            for (String tag : archive.get(UID).getTags()) {
                               if (tagList[i].equals(tag)) {
                                  tagCount++;
                               }
                            }
                         }
                    if (tagCount == tagList.length)
                                             return true;
                                        else
                                             return false;
                    }.collect()

Overlay overlay = new Overlay()
Random ran = new Random()

moleculesWithTags.stream().forEach{ UID ->
   Color color = new Color(ran.nextFloat(), ran.nextFloat(), ran.nextFloat())
   table = archive.get(UID).getDataTable()
   for (int i=0; i< table.getRowCount() ; i++ ) {
      PointRoi peakRoi = new PointRoi(table.getValue("x", i) + 0.5, table.getValue("y", i) + 0.5)
      //Have to select Use Names as Labels in Image>Overlay>Labels...
      peakRoi.setName(UID.substring(0, 5))
      peakRoi.setStrokeColor(color)
      //0-4 as options...
      peakRoi.setSize(4)
      //(0=hybrid, 1=crosshair, 2=dot, 3=circle)
      peakRoi.setPointType(2)
      peakRoi.setPosition((int)table.getValue("slice", i))
      overlay.add(peakRoi)
      }
}
image.setOverlay(overlay)
