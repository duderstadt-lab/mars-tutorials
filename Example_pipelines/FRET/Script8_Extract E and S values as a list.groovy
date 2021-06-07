#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

def E = []

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasParameter("E")){
  double Eval = molecule.getParameter("E")
  E.add(Eval)
  }
})

println(E)