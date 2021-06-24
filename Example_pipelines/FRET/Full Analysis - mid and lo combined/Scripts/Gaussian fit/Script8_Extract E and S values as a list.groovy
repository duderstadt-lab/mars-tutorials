#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

def E = []
def S = []

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasParameter("E")){
  double Eval = molecule.getParameter("E")
  double Sval = molecule.getParameter("S")
  E.add(Eval)
  S.add(Sval)
  }
})

println("E-values for this archive (FRET population)")
println(E)
println("")

println("S-values for this archive (FRET population)")
println(S)