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
//written by: prof. dr. Karl E. Duderstadt (Duderstadt lab)

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET


#@ Dataset acceptor_excitation_profile
#@ Dataset donor_excitation_profile
#@ MoleculeArchive archive
#@ OpService opService

import net.imglib2.view.Views

int maxAEX = opService.stats().max(acceptor_excitation_profile).getInteger()
int maxDEX = opService.stats().max(donor_excitation_profile).getInteger()

def aex_iMap_ra = Views.extendMirrorSingle(acceptor_excitation_profile.getImgPlus()).randomAccess()
def dex_iMap_ra = Views.extendMirrorSingle(donor_excitation_profile.getImgPlus()).randomAccess()

archive.molecules().forEach{ molecule ->
	molecule.getTable().rows().forEach{ row ->
		double aex = row.getValue("637_Red")
		int aex_x = (int) row.getValue("637_Red_X")
		int aex_y = (int) row.getValue("637_Red_Y")
		int dex_x = (int) row.getValue("532_Green_X")
		int dex_y = (int) row.getValue("532_Green_Y")

		double aex_corr = aex*(dex_iMap_ra.setPositionAndGet(dex_x, dex_y).getRealDouble()/maxDEX)/(aex_iMap_ra.setPositionAndGet(aex_x, aex_y).getRealDouble()/maxAEX)
		row.setValue("637_Red_Profile_Corrected", aex_corr)
	}
}
