/ $Id$
/ author:  JN Damask
/ descrip: Various tools to assist with the analysis of taq data. 

/ prints a logline
/ msg_: type string
.taq.logline: {[msg_]
  0N!(string .z.Z), "   taq |  ", msg_;
  };

/ returns bool. path_ is a string, e.g. "/home/user"
.taq.path_exists: {[path_]
  not () ~ key hsym "S"$ path_
  };

/ returns a bool. file_ is a string, e.g. "my_file.csv".
/   file_ is either in the current path or must be fully qualified:
/     "/home/user/data/my_file.csv"
.taq.file_exists: {[file_]
  not () ~ key hsym "S"$ file_
  };

/ saves a table as to a csv file. 
/ file_:  type string
/ table_: type table
.taq.save_csv: {[file_; table_]

  / left 0: right
  / right: .h.cd makes a comma-delimited string from the table
  / left: a file handle with name file_
  (hsym "S"$ file_) 0: .h.cd table_;

  };

/ import a taq quote csv file into kdb. 
/ file_ is a string.
.taq.import_quote_file: {[file_]

  if [not .taq.file_exists[file_];
    .taq.logline["file ", file_, " not found."];
    :()
  ];

  / load the quote file into a 'quote' table
  / the file must be formatted like:
  /  SYMBOL,DATE,TIME,BID,OFR,BIDSIZ,OFRSIZ,MODE,EX,MMID
  /  AA,20100105,9:30:00,16.76,16.88,4,1,12,Z,
  /  AA,20100105,9:30:00,16.81,16.84,6,2,12,T,
  /  AA,20100105,9:30:00,16.81,16.84,9,2,12,T,
  /  AA,20100105,9:30:00,16.81,16.84,12,2,12,T,
  /  ..
  `quote set
    ("SDTFFIIIC"; enlist ",") 0: hsym "S"$ file_;

  .taq.logline["loaded file ", file_];
  .taq.logline["  there are ", (string count quote), " records"];

  };

/ import a taq trade csv file into kdb. 
/ file_: type string.
.taq.import_trade_file: {[file_]

  if [not .taq.file_exists[file_];
    .taq.logline["file ", file_, " not found"];
    :()
  ];

  / load the trade file into a 'trade' table
  / the file must be formatted like:
  /  SYMBOL,DATE,EXCHANGE,TIME,PRICE,SIZE,COND,CORR,G127
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0
  /  CSCO,01/29/2010,P,7:39:34,22.6300,100,@F,0,0
  `trade set
    ("SDCTFIS"; enlist ",") 0: hsym "S"$ file_;

  .taq.logline["loaded file ", file_];
  .taq.logline["  there are ", (string count trade), " records"];

  };

/ makes a ruler in time (for one day) with intervals d_min minutes
/   apart. A table called 'ruler' is created.
/ start_: type time
/ end_:   type time
/ dmin_:  type int
.taq.make_time_ruler: {[start_; end_; dmin_]

  / convert to integers
  s_min: `int$ `minute$ start_;
  e_min: `int$ `minute$ end_;

  / find maximum number of intervals that fit the range
  n_intervals: ceiling (e_min - s_min) % dmin_;

  / intervals are marked from the end backwards to roughly
  / the start, and the start is explicitly added to the list. 
  time_v: `time$ `minute$ distinct s_min, reverse e_min - dmin_ * til n_intervals;

  / make a table called 'ruler' with column name TIME
  `ruler set 
    flip (enlist `TIME) ! enlist time_v;

  };

/ Given a quote table, a time ruler, symbol and exchange
/  returns a table of most-recent quotes as of the times
/  on the time ruler and adds the CNT column which is a 
/  count of the # of records between each time-point.
/ symbol_: type string
/ exch_: type string
/ time_ruler_: constructed from .taq.make_time_ruler[..]
.taq.make_quote_bars: {[symbol_; exch_; time_ruler_]

  / reorders the columns of the final table to that of quote
  ((cols quote), `CNT) xcols

    / joins time_ruler back into the result table
    / 'join' is the comma ,
    / 'each' is an adverb and is the single quote '
    / 'join-each' is ,'
    time_ruler_ ,'

      / take difference of CNT to get # quotes in each interval      
      update CNT: deltas CNT from 

        / asof join between selected quotes and the time_ruler
        /  the update adds the row index where the joins are
        /  made for the resulting table
        (update CNT:i from 
          select from quote where SYMBOL="S"$ symbol_, EX=exch_, MODE=12
        )
        asof
        time_ruler_
  };

/ Given a trade table, a time ruler, symbol and exchange
/  returns a table of most-recent trades as of the times
/  on the time ruler and adds the CNT and VOL columns.
/  CNT, like that for make_quote_bars, is the number of trade
/  events in each interval.
/  VOL is the accumulated traded volume per interval.
/ symbol_: type string
/ exch_: type string
/ time_ruler_: constructed from .taq.make_time_ruler[..]
.taq.make_trade_bars: {[symbol_; exch_; time_ruler_]

  / constrained selection from trade
  T: select from trade where SYMBOL="S"$ symbol_, EXCHANGE=exch_, COND in (`;`$"F";`$"@";`$"@F");

  / I don't need a local variable but the code is less cryptic this way.
  / creation of t follows the form in make_quote_bars[]
  t: ((cols trade), `CNT) xcols

    time_ruler_ ,'

      / note: unlike quotes I don't take CNT differences yet

      (update CNT:i from T) 
      asof
      time_ruler_;

  / eliminate null records if they exist
  t: delete from t where SYMBOL=`;
  
  / vector cut:
  /   list_l _ list_r
  / cuts list_r at indices specified by list_l, giving a list of lists
  / sum each (list of lists) sums the values of each list, giving a list
  update VOL: sum each t[`CNT] _ T[`SIZE],
         CNT: deltas CNT 
    from t
  };

/ taq timestamps are only resolved to the second.
/ this function takes the latest record per second (when it
/   exists) and appends the number of quote events in the
/   preceding second. 
/ symbol_: type string
/ exch_: type string
.taq.consolidate_quote_to_seconds: {[symbol_; exch_]

  / left ! right
  / right: a list of all seconds on the quote table
  / left: column name TIME 
  / x ! y is a dictionary, and flip x ! y is a table. 
  time_ruler: flip (enlist `TIME) ! (enlist asc exec distinct TIME from quote);

  .taq.make_quote_bars[symbol_; exch_; time_ruler]

  };

/ this function takes the latest record per second (when it
/   exists) and appends the number of trade events and 
/   traded volume within the second in the preceding second.
/ symbol_: type string
/ exch_: type string
.taq.consolidate_trade_to_seconds: {[symbol_; exch_]

  time_ruler: flip (enlist `TIME) ! (enlist asc exec distinct TIME from trade);

  .taq.make_trade_bars[symbol_; exch_; time_ruler]

  };
