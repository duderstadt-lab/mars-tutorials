//Generate a table with all values of a specified parameter 

#@ MoleculeArchive archive
#@ String parameterName
#@output MarsTable parameterTable

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

parameterTable = new MarsTable()
DoubleColumn col1 = new DoubleColumn(parameterName)

//Have to use this type for String columns...
GenericColumn col2 = new GenericColumn("UID")

parameterTable.add(col1)
parameterTable.add(col2)

archive.forEach{ molecule ->
    parameterTable.appendRow()
    parameterTable.setValue(parameterName, parameterTable.getRowCount() - 1, molecule.getParameter(parameterName))
    parameterTable.setValue("UID", parameterTable.getRowCount() - 1, molecule.getUID())
}
