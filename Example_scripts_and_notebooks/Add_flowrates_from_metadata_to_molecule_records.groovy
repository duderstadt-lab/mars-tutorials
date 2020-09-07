//Adds flowrates from metadata to molecule records

#@ MoleculeArchive archive

import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

archive.lock()

MarsMetadata metadata = archive.getMetadata(0)

archive.molecules().forEach{ molecule ->
  MarsTable table = molecule.getDataTable()
  if (!table.hasColumn("Flowrate (ul/min)"))
    table.appendColumn("Flowrate (ul/min)")

  for (int row=0;row<table.getRowCount();row++){
    table.set("Flowrate (ul/min)", row,metadata.getImage(0).getPlane(0, 0, row).getField("Flowrate (ul/min)") )
}}
---
// Add flowrates from a T,flowrate table to individual molecule records
#@ MarsTable tToFlowTable
#@ MoleculeArchive archive

import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*
import java.util.HashMap;

//First lets build a map from slice to flowrate
tToFlowMap = [:]

for (int row=0;row<tToFlowTable.getRowCount();row++)
   tToFlowMap.put((int) tToFlowTable.getValue("T", row), tToFlowTable.getValue("flowrate", row))

MarsMetadata metadata = archive.getMetadata(0)



archive.molecules().forEach{ molecule ->
  MarsTable table = molecule.getDataTable()
  if (!table.hasColumn("Flowrate (ul/min)"))
    table.appendColumn("Flowrate (ul/min)")

  for (int row=0;row<table.getRowCount();row++){
    table.set("Flowrate (ul/min)", row,metadata.getImage(0).getPlane(0, 0, row).getField("Flowrate (ul/min)") )
}}



archive.putMetadata(metadata)



---
//First make key value pairs of T values and flow rates
TToFlowrate = [:]


TToFlowrate.put((int) 1.0, 2.0)
TToFlowrate.put((int) 1.0, metadata.getImage(0).getPlane(0, 0, 1).getField("Flowrate (ul/min)"))

println(TToFlowrate)


for (int row=0;metadata.getImage(0).getSizeT();row++){
  TToFlowrate.put(row, metadata.getImage(0).getPlane(0, 0, row).getField("Flowrate (ul/min)") )
}

println(TToFlowrate)



println(metadata.getImage(0).getPlane(0, 0, 1).getField("Flowrate (ul/min)"))
println(metadata.getImage(0).getSizeT())


-------

sliceToFlowrate = [:]

for (int row=0;row<metaTable.getRowCount();row++) {
   sliceToFlowrate.put(metaTable.getValue("slice", row), metaTable.getValue("Flowrate (ul/min)", row))
}

println(sliceToFlowrate)

archive.getMoleculeUIDs().stream().forEach{ UID ->
   molecule = archive.get(UID)
   table = molecule.getDataTable()

   //If the column already exists we don't need to add it
   //instead we will just be overwriting the values below..
   if (!table.hasColumn("Flowrate (ul/min)"))
      table.appendColumn("Flowrate (ul/min)")

   for (int row=0;row<table.getRowCount();row++) {
      table.set("Flowrate (ul/min)", row, sliceToFlowrate.get(table.get("slice", row)))
   }

   archive.put(molecule)
}

archive.unlock()
