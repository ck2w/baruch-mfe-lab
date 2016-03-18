/ MTH9893 Time Series Spring 2016 Homework 1
/ Team: ShengQuan Zhou, Chi (Franklin) Ma

data_path: "/home/quan/Workspace/q_lab/data";
script_path: "/home/quan/Workspace/q_lab/scripts/q";

system "l ", script_path, "/taq_tools.q";

/ Problem #1: loading a TAQ quote file

.taq.logline["Problem #1:"];

/ a) Explain the code `quote set ("SDTFFIIIC"; enlist ",") 0: hsym "S"$ file_;
/    that appears in the function import_quote_file[].
/ Answer: The keyword 0: can be used in several situations:
/         1) Save text into file
/         2) Load fixed-width records
/         3) Load delimited records
/         See http://code.kx.com/wiki/Reference/ZeroColon for details.
/         Here 0: is used to load delimited records in the following way:
/         (types; delimiter) 0: filehandle
/         where in the function import_quote_file[]:
/               "SDTFFIIIC" is the type, enlist "," returns a list of delimiters,
/               "S"$ file_ converts the filename into a symbol,
/               and hsym converts its symbol argument into a file handle.
/         See http://code.kx.com/wiki/Reference/hsym for details.
/         Finally, the keyword set is an assignment operator, which assigns
/         the results loaded to the variable quote.

.taq.import_quote_file[data_path, "/quote/taq_20100105_quotes_dow30.csv"];

/ b) Number of records in table quote = 15089777
.taq.logline["Number of records in table quote:", string count quote];

/ c) What are the types of each column?

.taq.logline["Column types:"];
show meta quote

/ c     | t f a
/ ------| -----
/ SYMBOL| s    
/ DATE  | d    
/ TIME  | t    
/ BID   | f    
/ OFR   | f    
/ BIDSIZ| i    
/ OFRSIZ| i    
/ MODE  | i    
/ EX    | c  

/ d) Extract a list of the column names:

.taq.logline["Column names:", "," sv string cols quote];

/ Column names:SYMBOL,DATE,TIME,BID,OFR,BIDSIZ,OFRSIZ,MODE,EX"

/ e) Extract a list of unique SYMBOLS:

.taq.logline["A list of unique SYMBOLS:"];
distinct quote[`SYMBOL]

/ A list of unique SYMBOLS:
/ `AA`AXP`BA`BAC`CAT`CSCO`CVX`DD`DIS`GE`HD`HPQ`IBM`INTC`JNJ`JPM`KFT`KO`MCD`MMM`MRK`MSFT`PFE`PG`T`TRV`UTX`VZ`WMT`XOM

/ f) Query for the number of records for each SYMBOL:

.taq.logline["Number of records for each SYMBOL:"];
show select cnt:count i by SYMBOL from quote

/ Number of records for each SYMBOL:
/ SYMBOL| cnt    
/ ------| -------
/ AA    | 767141 
/ AXP   | 472911 
/ BA    | 283397 
/ BAC   | 1237552
/ CAT   | 358487 
/ CSCO  | 741070 
/ CVX   | 524716 
/ DD    | 353140 
/ DIS   | 395584 
/ GE    | 666941 
/ ..

/ g) Query for the number of exchanges (EX) for the SYMBOL AA:

show select cnt:count i by EX from quote where SYMBOL=`AA

/ EX| cnt   
/ --| ------
/ B | 78932 
/ C | 16251 
/ D | 9980  
/ I | 74164 
/ M | 5150  
/ N | 213820
/ P | 96639 
/ T | 175778
/ W | 8036  
/ Z | 88391 

/ h) Query for the number of modes (MODE) for the SYMBOL AA:

show select cnt:count i by MODE from quote where SYMBOL=`AA

/ MODE| cnt   
/ ----| ------
/ 5   | 3     
/ 10  | 1     
/ 12  | 767134
/ 15  | 3 

/ i) Find on the internet the association of these taq exchange
/    characters to the full name of the exchange. Also find, if
/    you can, the explanation of the MODE column.

/ Answer: http://www.nyxdata.com/data-products/daily-taq#155
/ A – NYSE MKT Stock Exchange
/ B – NASDAQ OMX BX Stock Exchange
/ C – National Stock Exchange
/ D – FINRA
/ I – International Securities Exchange
/ J – Direct Edge A Stock Exchange
/ K – Direct Edge X Stock Exchange
/ M – Chicago Stock Exchange
/ N – New York Stock Exchange
/ T – NASDAQ OMX Stock Exchange
/ P – NYSE Arca SM
/ S – Consolidated Tape System
/ T/Q – NASDAQ Stock Exchange
/ W – CBOE Stock Exchange
/ X – NASDAQ OMX PSX Stock Exchange
/ Y – BATS Y-Exchange
/ Z – BATS Exchange

/ The MODE column shows you what trading mode you are in.
/ 0: invalid field
/ 1: Slow Quote on the Offer
/ 2: Slow Quote on the Bid
/ 3: Closing Quote
/ 4: News dissemination
/ 5: Slow quote on the offer due to an LRP or Gap Quote
/ ..

/ Problem #2: loading a TAQ trade file

.taq.logline["Problem #2:"];

.taq.import_trade_file[data_path, "/trade/taq_ALL_20100105_trades.csv"];

/ a) What are the types of each column?

show meta trade

/ c       | t f a
/ --------| -----
/ SYMBOL  | s    
/ DATE    | d    
/ EXCHANGE| c    
/ TIME    | t    
/ PRICE   | f    
/ SIZE    | i    
/ COND    | s  


/ b) Extract a list of the column names:

.taq.logline["Column names:", "," sv string cols trade];

/ Column names:SYMBOL,DATE,EXCHANGE,TIME,PRICE,SIZE,COND

/ c) Query for the number of records for each SYMBOL:

show select cnt:count i by SYMBOL from trade

/ SYMBOL| cnt   
/ ------| ------
/ AA    | 119346
/ AXP   | 40162 
/ BA    | 38902 
/ BAC   | 270172
/ CAT   | 30093 
/ CSCO  | 118171
/ CVX   | 45046 
/ DD    | 38764 
/ DIS   | 39348 
/ GE    | 104651
/ ..

/ d) Query for the number of records for each SYMBOL, EXCHANGE:

show select cnt:count i by SYMBOL,EXCHANGE from trade

/ SYMBOL EXCHANGE| cnt  
/ ---------------| -----
/ AA     B       | 12751
/ AA     C       | 282  
/ AA     D       | 43094
/ AA     I       | 1696 
/ AA     M       | 44   
/ AA     N       | 10696
/ AA     P       | 15637
/ AA     T       | 21321
/ AA     W       | 36   
/ AA     Z       | 13789
/ ..

/ e) Find on the internet the association of the COND field with
/    condition description.

/ Answer: http://www.nyxdata.com/data-products/daily-taq#155
/ A = Slow on the Ask Side
/ B = Slow on the Bid Side
/ C = Closing
/ D = News Dissemination
/ E = Slow on the Bid due to LRP or GAP Quote
/ F = Slow on the Ask due to LRP or GAP Quote
/ G = Trading Range Indication
/ H = Slow on the Bid and Ask side
/ I = Order Imbalance
/ J = Due to a Related Security - News Dissemination
/ K = Due to a Related Security - News Pending
/ L = Closed Market Maker (NASD)
/ M = Additional Information
/ N = Non-firm quote
/ O = Opening Quote
/ P = News Pending
/ Q = Additional Information - Due to Related Security
/ R = Regular, two-sided open quote
/ S = Due to Related Security
/ T = Resume
/ U = Slow on the Bid and Ask due to LRP or GAP Quote
/ V = In View of Common
/ W = Slow Quote due to a Set Slow list on both the bid and offer sides
/ X = Equipment Changeover
/ Y = Regular - One Sided Quote (NASDAQ)
/ Z = No open/no resume  


/ Problem #3: Uniform Wall-Clock Time Ruler: Using the function
/             make_time_ruler[] as a template, write the function
/             make_time_sec_ruler[] to generate time points uniformly
/             spaced in seconds. The argument list of this function
/             should read: make_time_sec_ruler[start_; end_; dsec_].
/             Write this function into the .taq context.

.taq.logline["Problem #3:"];
.taq.make_time_sec_ruler: {[start_; end_; dsec_]

  / convert to integers
  s_sec: `int$ `second$ start_;
  e_sec: `int$ `second$ end_;

  / find maximum number of intervals that fit the range
  n_intervals: ceiling (e_sec - s_sec) % dsec_;

  / intervals are marked from the end backwards to roughly
  / the start, and the start is explicitly added to the list. 
  time_v: `time$ `second$ distinct s_sec, reverse e_sec - dsec_ * til n_intervals;

  / make a table called 'ruler' with column name TIME
  `ruler set 
    flip (enlist `TIME) ! enlist time_v;

  };


/ Problem #4: Prevaling Quote and Intensity Bins.
/             This problem is associated with the function make_quote_bars.

.taq.logline["Problem #4:"];

/ The data is all for the same date. For the time range, it is chosen to be 
/ 9:30-16:00, suggested by Professor Damask on the forum. 
/ The time interval is chosen to be 5 minutes.

.taq.make_time_ruler[09:30:00; 16:00:00; 5];
symbol_: "AA";
exch_: "T";
time_ruler_: ruler;

/ a) Execute the following query, report and explain your result
show update CNT:i from select from quote where SYMBOL="S"$ symbol_, EX=exch_, MODE=12

/ the row index i is appended to the results returned by select

/ SYMBOL DATE       TIME         BID   OFR   BIDSIZ OFRSIZ MODE EX CNT
/ --------------------------------------------------------------------
/ AA     2010.01.05 09:30:00.000 16.81 16.84 6      2      12   T  0  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 9      2      12   T  1  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 12     2      12   T  2  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 18     2      12   T  3  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 18     8      12   T  4  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 21     8      12   T  5  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 21     7      12   T  6  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 22     7      12   T  7  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 22     8      12   T  8  
/ AA     2010.01.05 09:30:00.000 16.81 16.84 23     8      12   T  9  
/ ..

/ b) Execute the following query, report and explain your result

t1: (update CNT:i from select from quote where SYMBOL="S"$ symbol_, EX=exch_, MODE=12)
    asof time_ruler_;

show t1

/ SYMBOL DATE       BID   OFR   BIDSIZ OFRSIZ MODE EX CNT  
/ ---------------------------------------------------------
/ AA     2010.01.05 16.82 16.84 15     6      12   T  20   
/ AA     2010.01.05 16.65 16.66 38     9      12   T  3835 
/ AA     2010.01.05 16.59 16.6  63     2      12   T  8667 
/ AA     2010.01.05 16.62 16.64 8      91     12   T  12845
/ AA     2010.01.05 16.63 16.64 61     103    12   T  16619
/ AA     2010.01.05 16.71 16.72 71     31     12   T  19302
/ AA     2010.01.05 16.6  16.61 6      91     12   T  23464
/ AA     2010.01.05 16.62 16.63 81     7      12   T  29672
/ AA     2010.01.05 16.57 16.58 145    78     12   T  31475
/ AA     2010.01.05 16.62 16.63 124    42     12   T  33674
/ ..

/ The latest-nearest record returned by each tick on the ruler.
/ Notice that the TIME column disappears because this is the column
/ being joined on by the time ruler.

/ Again, the cumulative count is appended to the right of the results
/ returned by the select. Because of the existence of the ruler,
/ the number of records between a successive pair of ruler ticks is
/ greater than 1.

/ c) Execute the query, report and explain your result

t2: update CNT: deltas CNT from t1

show t2

/ SYMBOL DATE       BID   OFR   BIDSIZ OFRSIZ MODE EX CNT 
/ --------------------------------------------------------
/ AA     2010.01.05 16.82 16.84 15     6      12   T  20  
/ AA     2010.01.05 16.65 16.66 38     9      12   T  3815
/ AA     2010.01.05 16.59 16.6  63     2      12   T  4832
/ AA     2010.01.05 16.62 16.64 8      91     12   T  4178
/ AA     2010.01.05 16.63 16.64 61     103    12   T  3774
/ AA     2010.01.05 16.71 16.72 71     31     12   T  2683
/ AA     2010.01.05 16.6  16.61 6      91     12   T  4162
/ AA     2010.01.05 16.62 16.63 81     7      12   T  6208
/ AA     2010.01.05 16.57 16.58 145    78     12   T  1803
/ AA     2010.01.05 16.62 16.63 124    42     12   T  2199
/ ..

/ deltas returns the successive difference between the CNT
/ column and the result is appended to the right of t1.

/ d) Execute the query, report and explain your result

show time_ruler_ ,' t2;

/ TIME         SYMBOL DATE       BID   OFR   BIDSIZ OFRSIZ MODE EX CNT 
/ ---------------------------------------------------------------------
/ 09:30:00.000 AA     2010.01.05 16.82 16.84 15     6      12   T  20  
/ 09:35:00.000 AA     2010.01.05 16.65 16.66 38     9      12   T  3815
/ 09:40:00.000 AA     2010.01.05 16.59 16.6  63     2      12   T  4832
/ 09:45:00.000 AA     2010.01.05 16.62 16.64 8      91     12   T  4178
/ 09:50:00.000 AA     2010.01.05 16.63 16.64 61     103    12   T  3774
/ 09:55:00.000 AA     2010.01.05 16.71 16.72 71     31     12   T  2683
/ 10:00:00.000 AA     2010.01.05 16.6  16.61 6      91     12   T  4162
/ 10:05:00.000 AA     2010.01.05 16.62 16.63 81     7      12   T  6208
/ 10:10:00.000 AA     2010.01.05 16.57 16.58 145    78     12   T  1803
/ 10:15:00.000 AA     2010.01.05 16.62 16.63 124    42     12   T  2199
/ ..

/ the asof join removes the column (TIME) being joined on, 
/ here the TIME column from the ruler is added back. 


/ e) Execute the command, report and explain your result

(cols quote), `CNT

/ `SYMBOL`DATE`TIME`BID`OFR`BIDSIZ`OFRSIZ`MODE`EX`CNT

/ (cols quote) returns the list of columns of table quote
/ (cols quote), `CNT appends `CNT to the end of the column list. 
/ The result is a list of column names with CNT appended at the end.

/ f) The function .taq.make_quote_bars[...] returns a table that is the result
/    of the preceding sequence of commands. Save an example of this table
/    to a csv file using .taq.save_csv["filename"; .taq.make_quote_bars[...]]
/    and plot CNT .vs. TIME. Add this plot to your homework.

show ((cols quote), `CNT) xcols time_ruler_ ,' t2

/ SYMBOL DATE       TIME         BID   OFR   BIDSIZ OFRSIZ MODE EX CNT 
/ ---------------------------------------------------------------------
/ AA     2010.01.05 09:30:00.000 16.82 16.84 15     6      12   T  20  
/ AA     2010.01.05 09:35:00.000 16.65 16.66 38     9      12   T  3815
/ AA     2010.01.05 09:40:00.000 16.59 16.6  63     2      12   T  4832
/ AA     2010.01.05 09:45:00.000 16.62 16.64 8      91     12   T  4178
/ AA     2010.01.05 09:50:00.000 16.63 16.64 61     103    12   T  3774
/ AA     2010.01.05 09:55:00.000 16.71 16.72 71     31     12   T  2683
/ AA     2010.01.05 10:00:00.000 16.6  16.61 6      91     12   T  4162
/ AA     2010.01.05 10:05:00.000 16.62 16.63 81     7      12   T  6208
/ AA     2010.01.05 10:10:00.000 16.57 16.58 145    78     12   T  1803
/ AA     2010.01.05 10:15:00.000 16.62 16.63 124    42     12   T  2199
/ ..

.taq.save_csv["quote_count.csv"; .taq.make_quote_bars[symbol_; exch_; ruler]];


/ Problem #5: Trade-Intensity Bins.
/             This problem is associated with the function make_trade_bars.

.taq.logline["problem #5:"];

.taq.make_time_ruler[09:30:00; 16:00:00; 5];

symbol_: "AA";
exch_: "T";
time_ruler_: ruler;

/ a) Execute the following query, report and explain your result

show update CNT:i from select from trade where SYMBOL="S"$ symbol_, EXCHANGE=exch_, COND in (`;`$"F";`$"@";`$"@F");

/ SYMBOL DATE       EXCHANGE TIME         PRICE SIZE COND CNT
/ -----------------------------------------------------------
/ AA     2010.01.05 T        07:41:57.000 16.6  500  F    0  
/ AA     2010.01.05 T        08:17:13.000 16.67 2480 F    1  
/ AA     2010.01.05 T        08:17:38.000 16.67 300  F    2  
/ AA     2010.01.05 T        08:19:45.000 16.67 498  F    3  
/ AA     2010.01.05 T        08:27:20.000 16.73 600  F    4  
/ AA     2010.01.05 T        08:30:42.000 16.74 1000 F    5  
/ AA     2010.01.05 T        08:30:42.000 16.74 1200 F    6  
/ AA     2010.01.05 T        08:30:47.000 16.74 500  F    7  
/ AA     2010.01.05 T        08:30:55.000 16.76 198  F    8  
/ AA     2010.01.05 T        08:30:55.000 16.76 400  F    9  
/ ..

/ All trading records for American Airline on NASDAQ with condition code F.
/ Append the row index to the right of the results returned by select.


/ b) Execute the following queries, report and explain your result

T: select from trade where SYMBOL="S"$ symbol_, EXCHANGE=exch_, COND in (`;`$"F";`$"@";`$"@F");
show T

/ SYMBOL DATE       EXCHANGE TIME         PRICE SIZE COND
/ -------------------------------------------------------
/ AA     2010.01.05 T        07:41:57.000 16.6  500  F   
/ AA     2010.01.05 T        08:17:13.000 16.67 2480 F   
/ AA     2010.01.05 T        08:17:38.000 16.67 300  F   
/ AA     2010.01.05 T        08:19:45.000 16.67 498  F   
/ AA     2010.01.05 T        08:27:20.000 16.73 600  F   
/ AA     2010.01.05 T        08:30:42.000 16.74 1000 F   
/ AA     2010.01.05 T        08:30:42.000 16.74 1200 F   
/ AA     2010.01.05 T        08:30:47.000 16.74 500  F   
/ AA     2010.01.05 T        08:30:55.000 16.76 198  F   
/ AA     2010.01.05 T        08:30:55.000 16.76 400  F   
/ ..

t: ((cols trade), `CNT) xcols time_ruler_ ,' (update CNT:i from T) asof time_ruler_;
show t

/ SYMBOL DATE       EXCHANGE TIME         PRICE SIZE COND CNT 
/ ------------------------------------------------------------
/ AA     2010.01.05 T        09:30:00.000 16.82 200  F    66  
/ AA     2010.01.05 T        09:35:00.000 16.66 100  F    642 
/ AA     2010.01.05 T        09:40:00.000 16.6  200  @    1248
/ AA     2010.01.05 T        09:45:00.000 16.63 100  @    1662
/ AA     2010.01.05 T        09:50:00.000 16.64 100  F    1966
/ AA     2010.01.05 T        09:55:00.000 16.71 800  F    2369
/ AA     2010.01.05 T        10:00:00.000 16.6  100  @    2867
/ AA     2010.01.05 T        10:05:00.000 16.63 100  @    3861
/ AA     2010.01.05 T        10:10:00.000 16.57 100  @    3969
/ AA     2010.01.05 T        10:15:00.000 16.63 500  F    4133
/ ..

/ ((cols trade), `CNT) is a list of columns of table trade with CNT appened.
/ Keyword xcols reorders columns in a table and returns a list of column names.
/ (update CNT:i from T) asof time_ruler_
/ returns the an additional CNT column which holds the cumulative count
/ as of the time ruler.

/ c) Execute the following commands, report and explain your result

t[`CNT]

/ 66 642 1248 1662 1966 2369 2867 3861 3969 4133 4404 4709 5130 5680 5888 6249 6784 7030 7355 7856 8518 8780 9234 9515 9746 9995 10399 10624 10964 11..
/ This command returns the CNT column of table t.

count t[`CNT]

/ 79
/ This command returns the number of records in column CNT, 
/ which equals the number of rows in table t.

T[`SIZE]

/ 500 2480 300 498 600 1000 1200 500 198 400 100 300 500 400 500 1000 100 100 400 500 100 100 500 500 100 400 100 100 100 400 102 4400 100 770 100 10..
/ This command returns the SIZE column of table T.

count T[`SIZE]

/ 21256
/ This command returns the number of records in column SIZE,
/ which equals the number of rows in table T.

t[`CNT] _ T[`SIZE]

/ 200 15656 300 100 100 100 350 250 200 100 250 100 100 100 100 100 100 100 100 100 100 100 600 100 100 100 200 300 300 160 300 100 100 100 300 200 ..
/ the keyword '_' represents a vector cut. 
/ Each entry in t[`CNT] is the index position in T[`SIZE] to cut and split
/ This option returns a splitted list of sublists of T[`SIZE],
/ each sublist containing the elements of trading volumes between cuts

sum each t[`CNT] _ T[`SIZE]

/ 123030 160939 112622 79532 86014 122049 260924 28142 48154 84727 74071 123810 219380 47690 89104 184180 65584 95512 165078 196063 95614 156265 8546..
/ returns the summation of trade volume (sum over SIZE*CNT) between each cut, 
/ which is a vector t[`CNT] long of accumulated trade sizes (or volumes). 


/ d) Execute the following query, report and explain your result

show update VOL: sum each t[`CNT] _ T[`SIZE], CNT: deltas CNT from t

/ SYMBOL DATE       EXCHANGE TIME         PRICE SIZE COND CNT VOL   
/ ------------------------------------------------------------------
/ AA     2010.01.05 T        09:30:00.000 16.82 200  F    66  123030
/ AA     2010.01.05 T        09:35:00.000 16.66 100  F    576 160939
/ AA     2010.01.05 T        09:40:00.000 16.6  200  @    606 112622
/ AA     2010.01.05 T        09:45:00.000 16.63 100  @    414 79532 
/ AA     2010.01.05 T        09:50:00.000 16.64 100  F    304 86014 
/ AA     2010.01.05 T        09:55:00.000 16.71 800  F    403 122049
/ AA     2010.01.05 T        10:00:00.000 16.6  100  @    498 260924
/ AA     2010.01.05 T        10:05:00.000 16.63 100  @    994 28142 
/ AA     2010.01.05 T        10:10:00.000 16.57 100  @    108 48154 
/ AA     2010.01.05 T        10:15:00.000 16.63 500  F    164 84727 
/ ..

/ The deltas CNT was executed in the make_quote_bars function, but here it is executed in this last update.

/ e) The function .taq.make_trade_bars[...] returns a table that is the result
/    of the preceding sequency of commands. Save an example of this table to a
/    csv file using .taq.save_csv["filename"; .taq.make_trade_bars[...]] and
/    plot VOL, CNT .vs. TIME. Add this plot to your homework.


.taq.save_csv["trade_volume.csv"; .taq.make_trade_bars[symbol_; exch_; ruler]];

