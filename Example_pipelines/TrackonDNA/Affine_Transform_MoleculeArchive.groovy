#@ MoleculeArchive archive

import net.imglib2.realtransform.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*

AffineTransform2D xfm = new AffineTransform2D()
xfm.set(1.003149836793469, 6.13768352145E-4, 1.223492388374801, 3.61289752847E-4, 1.002457751873398, 507.24885432325897)

archive.getMoleculeUIDs().stream().forEach{ UID ->
	Molecule molecule = archive.get(UID)
	MarsTable table = molecule.getTable()

	for (int row=0; row <table.getRowCount() ; row++ ) {
		double[] source = new double[2]
		source[0] = table.getValue("x", row)
		source[1] = table.getValue("y", row)

		double[] target = new double[2]

		xfm.apply(source, target)

		table.setValue("x", row, target[0])
		table.setValue("y", row, target[1])
	}

	archive.put(molecule)
}

archive.logln(LogBuilder.buildTitleBlock("Transform Molecule Coordinates"))
archive.logln("[m00, m01, m02, m10, m11, m12]")
archive.logln("[" +
xfm.get(0,0) + ", " + xfm.get(0,1) + ", " + xfm.get(0,2) + ", " +
xfm.get(1,0) + ", " + xfm.get(1,1) + ", " + xfm.get(1,2) + "]")
archive.logln(LogBuilder.endBlock())
