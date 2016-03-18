/ $Id$
/ author:  JN Damask
/ descrip: Examples using the taq_tools.q functions
/ use:     start q using either
/            $ rlwrap q
/          or (pref)
/            $ rlwrap q -p 18001
/          alter this script file for the quote date you want
/          load this script from the q prompt
/            q)\l taq_quote_examples.q
/          make your own scripts using these functions, as you like.

/ NOTE: I have converted the dow-30 quote files to one file per day,
/       each day file having all 30 stocks.

/ specify date and root path
taq_date: "20100105";
taq_path: "/home/quan/Workspace/q_lab";
taq_bar: 1;

/ import the tools script -- must specify path
@[system; "l ", taq_path, "/scripts/q/taq_tools.q"; {0N!"no good"; exit -1}];

/ import a quote file -- must specify path
/ saves import file to the 'quote' table
.taq.logline["loading quote data"];
.taq.import_quote_file[taq_path, "/data/quote/taq_", taq_date, "_quotes_dow30.csv"];

/ import a trade file -- must specify path
.taq.logline["loading trade data"];
.taq.import_trade_file[taq_path, "/data/trade/taq_ALL_", taq_date, "_trades.csv"];

/ make a time ruler for the bars
/ save the ruler to the 'ruler' table
.taq.make_time_ruler[09:30:00; 16:00:00; taq_bar];

/ iterate over all tickers in quote, for each ticker sample the quote table
/ then unlist (raze) the result into one table calls 'quote_bars'
.taq.logline["making quote bars on ", (string taq_bar), " min intervals"];
quote_bars: 
  raze
    {[s]
      .taq.make_quote_bars[string s; "T"; ruler]
    } each exec distinct SYMBOL from quote;

.taq.logline["  there are ", (string count quote_bars), " records in quote_bars"];

/ save the result to a csv file -- must specify path
.taq.fn: taq_path, "/data/taq_", taq_date, "_quote_", (string taq_bar), "_bars_dow30.csv";
.taq.logline["saving file ", .taq.fn];
.taq.save_csv[.taq.fn; quote_bars];

/ same for trades
.taq.logline["making trade bars on ", (string taq_bar), " min intervals"];
trade_bars:
  raze
    {[s]
      .taq.make_trade_bars[string s; "T"; ruler]
    } each exec distinct SYMBOL from trade;
.taq.logline["  there are ", (string count trade_bars), " records in trade_bars"];

/ save the result to a csv file -- must specify path
.taq.fn: taq_path, "/data/taq_", taq_date, "_trade_", (string taq_bar), "_bars_dow30.csv";
.taq.logline["saving file ", .taq.fn];
.taq.save_csv[.taq.fn; trade_bars];

