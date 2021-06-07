#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","0","active")){
  MarsTable table = molecule.getSegmentsTable("T","0","active")
  double ymax = table.max("A")
  molecule.setParameter("Iaemaex",ymax)}
  if (molecule.hasSegmentsTable("T","1 Green","active")){
  MarsTable table2 = molecule.getSegmentsTable("T","1 Green","active")
  double ymax2 = table2.max("A")
  molecule.setParameter("Idemdex",ymax2)}
  if (molecule.hasSegmentsTable("T","1 Red","active")){
  MarsTable table3 = molecule.getSegmentsTable("T","1 Red","active")
  double ymax3 = table3.max("A")
  molecule.setParameter("Iaemdex",ymax3)}
  }) 
