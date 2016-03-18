/ $Id$
/ author:  JN Damask
/ descrip: Example code for lecture 1

/ preamble
\c 15 150
data_path:   "/home/jaydamask/vm_share/teaching/Baruch/time_series/data";
script_path: "/home/jaydamask/vm_share/teaching/Baruch/time_series/scripts/q";
system "l ", script_path, "/taq_tools.q";

/----------------------------------------------------------
/ -- Loading a flat table

  /  SYMBOL,DATE,EXCHANGE,TIME,PRICE,SIZE,COND,CORR,G127
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0

  /  `trade set
  /    ("SDCTFIS"; enlist ",") 0: hsym "S"$ file_;

.taq.import_trade_file[data_path, "/trade/taq_ALL_20100105_trades.csv"];

show trade
show reverse trade

/----------------------------------------------------------
/ -- Flat table inspection

/ trade column names
cols trade

/ trade column types
meta trade

/ # of records
count trade

/ implicit index
update seq:i from trade

/ flip trade to see dict
flip trade

/ extract a column-value list via dictionary syntax
(flip trade)[`DATE]
trade[`DATE]
trade[`a]  <-  `a not in key set, returns null

/ extract the column names -- keys -- from a flat table
cols trade
trade[cols trade]

/----------------------------------------------------------
/ -- Make a table directly from a dictionary

D: () ! ();
D[`event]: (0;1;2);
D[`px]: (80.1; 80.2; 80.3);
D[`sz]: (100; 200; 100);

show D
show flip D
T: flip D;

meta T

add (k,v) pair that is jagged
D[`valid]: (1b; 1b; 1b; 0b);
flip D

/ correct this
D[`valid]: (1b; 1b; 0b);
flip D


/----------------------------------------------------------
/ -- Making distinct lists (sets)

symbol_set: distinct trade[`SYMBOL]
exch_set: distinct trade[`EXCHANGE]
cond_set: distinct trade[`COND]

/----------------------------------------------------------
/ -- Notes about lists

/ make a list having same type
L1: `ibm`csco`dtv`appl`goog

/ 'join' command, binary operator
L1 , `qqq

/ join two lists
L2: `xle`xlb
L1 , L2

/ note: a list of lists:
(L1; L2)

/ 'raze' command, unary opeator
raze (L1; L2)

/ 'cut' command, binary operator
L3: til 10
3 _ L3

/ cut is a vector operator, so:
(0;3) _ L3
(0;2;6;9) _ L3
/ (0;5;15) _ L3  / this gives an error


/----------------------------------------------------------
/ -- Select query

/ 'select c from T' form
select SYMBOL, DATE from trade
count select SYMBOL, DATE from trade

select CNT:i from trade

select min TIME, max TIME from trade

/ constraints
select from trade where SYMBOL=`AA
count select from trade where SYMBOL=`AA

/ 'where' keyword, unary operator
where trade[`EXCHANGE] = "P"

/ 'select c from T where k'
select from trade where SYMBOL=`AA, EXCHANGE = "P"
count select from trade where SYMBOL=`AA, EXCHANGE = "P"

/ careful about implicit index
select from (update CNT:i from trade) where SYMBOL=`AA, EXCHANGE="P"
update CNT:i from select from trade where SYMBOL=`AA, EXCHANGE="P"


/----------------------------------------------------------
/ -- Primary-Keyed Tables

/ dictionary of two flat tables, the key table having unique entries
kt: flip (enlist `EXCHANGE) ! (enlist `P`T);
vt: flip (`stocks_listed`comp_vol) ! ((2500;3400); (15;23));

show kt
show vt

show kt ! vt

ET: kt ! vt

/ key and value of p-key table
show key ET
show value ET

/ upsert
rec_1: (`Q; 3400; 9);
rec_2: (`P; 1200; 14);

show rec_1
show rec_2

show ET upsert rec_1
show ET upsert rec_2

/ a p-key table is itself NOT a dictionary in the sense of a flat table
/ flip ET

/ shorthand to remove a key from a table
show 0 ! ET
(0!ET)[`EXCHANGE]

/----------------------------------------------------------
/ -- Aggregation with select

/ 'select c by g from T'
select cnt: count i by EXCHANGE from trade

/ compare
select cnt: count i by SYMBOL, EXCHANGE from trade
select cnt: count i by EXCHANGE from trade where SYMBOL=`AA


/----------------------------------------------------------
/ -- Asof joins

AT: flip (`mark`val) ! ((til 20); sums rand each 20 # 1f);
M: flip (enlist `mark) ! enlist (3; 6; 10; 15; 22);

show AT
show M

show AT asof M
show M ,' AT asof M

/----------------------------------------------------------
/ -- Wrap up

/ casts $
`int $ 3.1415
`int $ 09:30:00.000
`second $ 09:30:00.000

`$"a_symbol"
"S"$"ya_symbol"

/ functions -- lambda calculus
f: {[v] 2 * v }
g: {[v] 3 * f[v] }

/ function notes: 
/  a) multi-line functions cannot be cut and paste directly
/  b) multi-line functions must have the closing '}' whitespaced away from 
/        the starting column
/  c) each line in a function must end w/ and ';'
/  d) if the last line ends in ';' then no return,
/     else the value of the last line is returned











