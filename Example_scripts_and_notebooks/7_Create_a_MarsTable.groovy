#@output MarsTable table

import de.mpg.biochem.mars.table.*;

//Initialize a new ResultsTable with 3 columns and no rows
table = new MarsTable(3,0)
table.setColumnHeader(0, "column1")
table.setColumnHeader(1, "column2")
table.setColumnHeader(2, "column3")

for (int row=0;row<10;row++) {
   //Since we initialised the table with no rows
   //We have to increase the table size by one row
   //before we can add more values
   table.appendRow()
   table.setValue("column1", row, row*0.25)
   table.setValue("column2", row, row*0.25)
   table.setValue("column3", row, row*0.25)
}

//Alternatively one can build the MarsTable in the following way

#@output MarsTable table

import de.mpg.biochem.mars.table.*
import org.scijava.table.*


//Initialize a new empty MARSResultsTable
table = new MarsTable()
DoubleColumn col1 = new DoubleColumn("column1")
DoubleColumn col2 = new DoubleColumn("column2")
DoubleColumn col3 = new DoubleColumn("column3")

for (int row=0;row<10;row++) {
   col1.add((double)row*0.25)
   col2.add((double)row*0.25)
   col3.add((double)row*0.25)
}

table.add(col1)
table.add(col2)
table.add(col3)
