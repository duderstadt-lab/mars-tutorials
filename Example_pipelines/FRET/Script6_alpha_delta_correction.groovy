/*******************************************************************************
 * Copyright (C) 2021, Duderstadt Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
//written by: Nadia M. Huisjes, MSc. (Duderstadt lab)

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET

// Calculate the alpha and delta correction factors and apply the correction on E and S

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

def E_list = []
def S_list = []

//calculate the corrected S and E values for the AO and DO populations respectively
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

//calculate the mean S and E values from the AO and DO populations and assign a parameter in the metadata
double E_mean = E_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/E_list.size()
double S_mean =  S_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/S_list.size()

archive.metadata().forEach{metadata -> 
	metadata.setParameter("iiEappDO_mean", E_mean)
	metadata.setParameter("iiSappAO_mean", S_mean)
}


//calculate alpha and delta and assign a parameter in the metadata
double alpha = E_mean/(1- E_mean)
double delta = S_mean/(1-S_mean)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("alpha", alpha)
	metadata.setParameter("delta", delta)
}

//calculate FAD for each FRET molecule and the corrected E and S values
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








