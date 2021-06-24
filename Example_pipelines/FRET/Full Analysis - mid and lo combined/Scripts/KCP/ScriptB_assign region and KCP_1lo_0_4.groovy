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

//Assign the global region "active"
#@ MoleculeArchive archive
#@ ImageJ ij

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

archive.metadata().forEach{metadata ->
	metadata.putRegion(new MarsRegion("active", "T", 2, 500, "#416EF4", 0.2))
	metadata.addTag("lo")
}

//KCP analysis of 1 Red vs. T
import de.mpg.biochem.mars.kcp.commands.*

final KCPCommand kcpCalc = new KCPCommand();

kcpCalc.setContext(ij.getContext());

kcpCalc.setArchive(archive);
kcpCalc.setXColumn("T");
kcpCalc.setYColumn("1 Red");
kcpCalc.setConfidenceLevel(0.99);
kcpCalc.setGlobalSigma(70000);
kcpCalc.setRegionSource("Metadata");
kcpCalc.setCalculateBackgroundSigma(false);
kcpCalc.setBackgroundRegion("None");
kcpCalc.setAnalyzeRegion(true);
kcpCalc.setRegion("active");
kcpCalc.setFitSteps(true);
kcpCalc.setIncludeTags("All");
kcpCalc.setTags("None");
kcpCalc.setThreads(8)

kcpCalc.run();

//KCP analysis of 0 vs. T
final KCPCommand kcpCalc2 = new KCPCommand();

kcpCalc2.setContext(ij.getContext());

kcpCalc2.setArchive(archive);
kcpCalc2.setXColumn("T");
kcpCalc2.setYColumn("0");
kcpCalc2.setConfidenceLevel(0.99);
kcpCalc2.setGlobalSigma(70000);
kcpCalc2.setRegionSource("Metadata");
kcpCalc2.setCalculateBackgroundSigma(false);
kcpCalc2.setBackgroundRegion("None");
kcpCalc2.setAnalyzeRegion(true);
kcpCalc2.setRegion("active");
kcpCalc2.setFitSteps(true);
kcpCalc2.setIncludeTags("All");
kcpCalc2.setTags("None");
kcpCalc2.setThreads(8)

kcpCalc2.run();

//KCP analyisis of 1 Green vs. T
final KCPCommand kcpCalc3 = new KCPCommand();

kcpCalc3.setContext(ij.getContext());

kcpCalc3.setArchive(archive);
kcpCalc3.setXColumn("T");
kcpCalc3.setYColumn("1 Green");
kcpCalc3.setConfidenceLevel(0.99);
kcpCalc3.setGlobalSigma(70000);
kcpCalc3.setRegionSource("Metadata");
kcpCalc3.setCalculateBackgroundSigma(false);
kcpCalc3.setBackgroundRegion("None");
kcpCalc3.setAnalyzeRegion(true);
kcpCalc3.setRegion("active");
kcpCalc3.setFitSteps(true);
kcpCalc3.setIncludeTags("All");
kcpCalc3.setTags("None");
kcpCalc3.setThreads(8)

kcpCalc3.run();