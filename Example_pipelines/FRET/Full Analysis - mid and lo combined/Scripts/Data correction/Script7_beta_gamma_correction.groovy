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

// Calculate the beta and gamma correction factors and apply the correction on E and S

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import org.apache.commons.math3.stat.regression.SimpleRegression

def regression = new SimpleRegression(true)
def E_list = []
def S_list = []

//create a list with all data points (iiiEapp and iiiSapp) to be used for the fit
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double E = molecule.getParameter("iiiEapp")
    double S = 1 / molecule.getParameter("iiiSapp")
    S_list.add(S)
    E_list.add(E)}})


//add the datapoints to be included in the regression analysis and store found parameters
for (i in [E_list,S_list].transpose()){
	regression.addData(i[0],i[1])
}

double a = regression.getSlope()
double b = regression.getIntercept()
double beta = a + b - 1
double gamma = (b - 1) / (a + b - 1)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("beta", beta)
	metadata.setParameter("gamma", gamma)
}

//calculate F_AA and F_DD values for each molecule and calculate the E and S values
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double FAD = molecule.getParameter("F_AD")
  	double Idemdex = molecule.getParameter("iiIdemdex")
  	double Iaemaex = molecule.getParameter("iiIaemaex")
  	double FDD = gamma * Idemdex
  	double FAA = Iaemaex / beta
    molecule.setParameter("F_DD",FDD)
    molecule.setParameter("F_AA",FAA)
    double E = FAD / (FDD + FAD)
    double S = (FDD + FAD) / (FDD + FAD + FAA)
    molecule.setParameter("E",E)
    molecule.setParameter("S",S)
  }
    
 })








    






