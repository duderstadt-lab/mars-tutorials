//Get the intensity of the first frame for all Molecules
//Change "green" to the name of the intensity column

#@ MoleculeArchive archive
#@output MarsTable table

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

table = new MarsTable()
DoubleColumn col1 = new DoubleColumn("Intensity")

//Have to use this type for String columns...
GenericColumn col2 = new GenericColumn("UID")

table.add(col1)
table.add(col2)

archive.forEach{ molecule ->
    table.appendRow()
    table.setValue("Intensity", table.getRowCount() - 1, molecule.getDataTable().getValue("Green", 0))
    table.setValue("UID", table.getRowCount() - 1, molecule.getUID())
} 
