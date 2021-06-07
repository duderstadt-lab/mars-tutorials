#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  double Iaemdex = molecule.getParameter("Iaemdex")
  double Idemdex = molecule.getParameter("Idemdex")
  double Iaemaex = molecule.getParameter("Iaemaex")
  double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  double E = Iaemdex / (Iaemdex + Idemdex)
  molecule.setParameter("iSapp",S)
  molecule.setParameter("iEapp",E)}
  })