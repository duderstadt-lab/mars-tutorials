#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

def E_list = []
def S_list = []

//Calculate the corrected S and E values for the AO and DO populations respectively
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("AO_active")){
  	double Iaemdex = molecule.getParameter("iiIaemdex")
    double Idemdex = molecule.getParameter("iiIdemdex")
    double Iaemaex = molecule.getParameter("iiIaemaex")
  	double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  	double E = (Iaemdex) / (Iaemdex + Idemdex)
    molecule.setParameter("iiSappAO",S)
    molecule.setParameter("iiEappAO",E)
    S_list.add(S)
  	}
  	else if (molecule.hasTag("DO_active")){
		double Iaemdex = molecule.getParameter("iiIaemdex")
        double Idemdex = molecule.getParameter("iiIdemdex")
        double Iaemaex = 0
  		double E = Iaemdex / (Iaemdex + Idemdex)
  		double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  		molecule.setParameter("iiEappDO",E)
  		molecule.setParameter("iiSappDO",S)
  		E_list.add(E)
  	}
  })

//Calculate the mean S and E values from the AO and DO populations and assign a parameter in the metadata
double E_mean = E_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/E_list.size()
double S_mean =  S_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/S_list.size()

archive.metadata().forEach{metadata -> 
	metadata.setParameter("iiEappDO_mean", E_mean)
	metadata.setParameter("iiSappAO_mean", S_mean)
}


//Calculate alpha and delta and assign a parameter in the metadata
double alpha = E_mean/(1- E_mean)
double delta = S_mean/(1-S_mean)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("alpha", alpha)
	metadata.setParameter("delta", delta)
}

//Calculate FAD for each FRET molecule and the corrected E and S values
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double Iaemdex = molecule.getParameter("iiIaemdex")
    double Idemdex = molecule.getParameter("iiIdemdex")
    double Iaemaex = molecule.getParameter("iiIaemaex")
    double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
    double iiiE = FAD / (FAD + Idemdex)
    double iiiS = (FAD + Idemdex) / (FAD + Idemdex + Iaemaex)
    molecule.setParameter("F_AD",FAD)
    molecule.setParameter("iiiEapp",iiiE)
    molecule.setParameter("iiiSapp",iiiS)
  }
})








