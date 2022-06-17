/*******************************************************************************
 * Copyright (C) 2022, Duderstadt Lab
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
//written by: Karl E. Duderstadt

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET

#@ String (label="Aem|Aex (format: channel_region)", value="637_Red") aemaex
#@ String (label="Dem|Dex (format: channel_region)", value="532_Green") demdex
#@ ImgPlus (label="Acceptor excitation profile") acceptor_excitation_profile
#@ ImgPlus (label="Donor excitation profile") donor_excitation_profile
#@ MoleculeArchive archive
#@ OpService opService

import net.imglib2.view.Views
import de.mpg.biochem.mars.util.*

//Build log message
builder = new LogBuilder()
String log = LogBuilder.buildTitleBlock("FRET workflow 2 profile correction")
builder.addParameter("Workflow version", "0.1")
builder.addParameter("Aem|Aex", aemaex)
builder.addParameter("Dem|Dex", demdex)
builder.addParameter("Acceptor excitation profile", acceptor_excitation_profile.getName())
builder.addParameter("Donor excitation profile", donor_excitation_profile.getName())
log += builder.buildParameterList()

int maxAEX = opService.stats().max(acceptor_excitation_profile).getInteger()
int maxDEX = opService.stats().max(donor_excitation_profile).getInteger()

def aex_iMap_ra = Views.extendMirrorSingle(acceptor_excitation_profile).randomAccess()
def dex_iMap_ra = Views.extendMirrorSingle(donor_excitation_profile).randomAccess()

archive.molecules().forEach{ molecule ->
	molecule.getTable().rows().forEach{ row ->
		double aex = row.getValue(aemaex)
		int aex_x = (int) row.getValue(aemaex + "_X")
		int aex_y = (int) row.getValue(aemaex + "_Y")
		int dex_x = (int) row.getValue(demdex + "_X")
		int dex_y = (int) row.getValue(demdex + "_Y")

		double aex_corr = aex*(dex_iMap_ra.setPositionAndGet(dex_x, dex_y).getRealDouble()/maxDEX)/(aex_iMap_ra.setPositionAndGet(aex_x, aex_y).getRealDouble()/maxAEX)
		row.setValue(aemaex + "_Profile_Corrected", aex_corr)
	}
}

log += "\n" + LogBuilder.endBlock(true) + "\n"
archive.logln(log)
