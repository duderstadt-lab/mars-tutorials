//Calculate the slope between the start and end of a trace
//Allows to set for start and end, x and y, and tag information
//Generates an output table and adds the calculated slope as molecule parameter

#@ MoleculeArchive archive
#@ String slopeColumn
#@ String xColumn
#@ double start
#@ double end
#@ String filterTag
#@output MarsTable tableOUT

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

tableOUT = new MarsTable()
DoubleColumn col1 = new DoubleColumn("slope")

GenericColumn col2 = new GenericColumn("UID")

tableOUT.add(col1)
tableOUT.add(col2)

archive.lock()

archive.getMoleculeUIDs().stream().filter({UID -> archive.get(UID).hasTag(filterTag)}).forEach({ UID ->
      Molecule molecule = archive.get(UID)
      MarsTable table = molecule.getDataTable()

     double[] fit = table.linearRegression(xColumn, slopeColumn, start, end)

     molecule.setParameter("slope", fit[2])

     tableOUT.appendRow()
     tableOUT.setValue("slope", tableOUT.getRowCount() - 1, fit[2])
     tableOUT.setValue("UID", tableOUT.getRowCount() - 1, UID)

      archive.put(molecule);
   })

archive.unlock()
