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

#@ String (label="Acceptor column (format: channel_region)", value="637_Red") acceptorColumn
#@ String (label="Donor column (format: channel_region)", value="532_Green") donorColumn
#@ MoleculeArchive archive
#@ ImageJ ij

import de.mpg.biochem.mars.kcp.commands.*
import de.mpg.biochem.mars.util.*

//Build log message
builder = new LogBuilder()
String log = LogBuilder.buildTitleBlock("FRET workflow 3 find bleaching positions")
builder.addParameter("Workflow version", "0.2")
builder.addParameter("Acceptor column", acceptorColumn)
builder.addParameter("Donor column", donorColumn)
log += builder.buildParameterList()
archive.logln(log)

//Make an instance of the Command you want to run
final SingleChangePointFinder scpCalc = new SingleChangePointFinder()

//Populates @Parameters Services etc. using the current context which we get from the ImageJ Input
scpCalc.setContext(ij.getContext())

//Set all the input parameters
scpCalc.setArchive(archive)
scpCalc.setXColumn("T")
scpCalc.setFitSteps(true)
scpCalc.setAddPosition(true)
scpCalc.setIncludeTags("Tagged with")
scpCalc.setThreads(8)

//Set to not use
scpCalc.setAnalyzeRegion(false)
scpCalc.setRegionSource("Molecules")
scpCalc.setRegion("")
scpCalc.setAddSegmentsTable(false)

//Add Acceptor_Bleach to FRET molecules
scpCalc.setPosition("Acceptor_Bleach")
scpCalc.setYColumn(acceptorColumn)
scpCalc.setTags("FRET")
scpCalc.run()

//Add Acceptor_Bleach to AO molecules
scpCalc.setTags("AO")
scpCalc.run()

//Add Donor_Bleach to FRET molecules
scpCalc.setPosition("Donor_Bleach")
scpCalc.setYColumn(donorColumn)
scpCalc.setTags("FRET")
scpCalc.run()

//Add Acceptor_Bleach to AO molecules
scpCalc.setTags("DO")
scpCalc.run()
