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

// Calculate the uncorrected apparent E and S values (iEapp and iSapp)

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

//retrieve the intensity parameters and calculate the iEapp and iSapp value for each molecule
//store the iEapp and iSapp values as a parameter in the respective molecule record
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