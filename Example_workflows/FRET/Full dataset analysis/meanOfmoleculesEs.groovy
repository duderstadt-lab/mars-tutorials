#@ MoleculeArchive archive

import org.scijava.table.*
import de.mpg.biochem.mars.table.*

DoubleColumn colElo = new DoubleColumn("1-lo E")
archive.molecules().filter{ m -> m.hasTag("FRET")}.filter{ m -> m.hasTag("Accepted")}.filter{ m -> m.hasTag("1-lo")}.forEach{ m ->
	colElo.add(m.getParameter("E"))
}
MarsTable tableLo = new MarsTable("1-lo table")
tableLo.add(colElo)
println "1-lo E " + tableLo.mean("1-lo E") + " std " + tableLo.std("1-lo E") + " sem " + tableLo.sem("1-lo E")



DoubleColumn colEmid = new DoubleColumn("1-mid E")
archive.molecules().filter{ m -> m.hasTag("FRET")}.filter{ m -> m.hasTag("Accepted")}.filter{ m -> m.hasTag("1-mid")}.forEach{ m ->
	colEmid.add(m.getParameter("E"))
}

MarsTable tableMid = new MarsTable("1-mid table")
tableMid.add(colEmid)
println "1-mid E " + tableMid.mean("1-mid E") + " std " + tableMid.std("1-mid E") + " sem " + tableMid.sem("1-mid E")