//Filter to tag all molecules where "reversal"<0.5 and tag those molecules as "background"

#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*

archive.stream().filter{ molecule -> molecule.getParameter("reversal") < 0.5}.forEach{ molecule ->
   //This would tag all molecule with reversal below 0.5
   molecule.addTag("background")
   archive.put(molecule);
}
