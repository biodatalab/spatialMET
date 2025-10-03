#!/usr/bin/perl -w

# 2025-03-25  speed up quote/space stripping
# 2025-03-22  slower, but more correct, malformed fields quote/space stripping
# 2025-03-21  bugfix --escape-dq
# 2025-03-21  fix enclosing and line wrap dequoting/despacing
# 2025-03-20  more support for unicode EOL-like and unicode spaces
# 2025-03-20  improve default unwrap spacing
# 2025-03-20  better handle potential presence/absence of [\r\n] in $ regex
# 2025-03-19  bugfix: line unwrapped returned string matches stdout output now
# 2025-03-19  add --escape-eol, improve default unwrap spacing
# 2025-03-18  --unwrap now handles lone closing " on line followed by EOL
# 2024-10-17  escape all $delim within regular expressions
# 2024-10-16  add --semicolon option for semicolon-delimited
# 2024-07-12  add Usage statement
# 2024-03-20  replace line wrapped EOL with spaces or empty
# 2024-03-19  add --escape-dq flag
# 2024-03-19  add support for line wraps with new --unwrap flag
# 2024-03-19  fix edge cases such as ^"a,b" text$ counting as enclosed field
# 2022-10-24  fix "commments" typo in 2022-07-12 changelog (oh, the irony)
# 2022-07-12  minor typo corrections in documentation comments
# 2021-12-09  minor typo corrections in documentation comments
# 2021-08-19  improve/update documentation of csv2tsv_not_excel() function
# 2021-08-16  strip leading/trailing spaces from enclosed fields,
#             instead of the previous intentional preservation behavior
# 2021-08-16  bugfix removal of enclosing double-quotes with trailing spaces
# 2021-06-21  special case "" in a field by itself
# 2020-08-05  modify split regex to handle a few more invalid edge cases
# 2020-08-04  modify split regex to handle a few more invalid edge cases
# 2020-07-31  more improvements to robustness and optimizations
# 2020-07-28  speed up by replacing more quick and easy cases first
# 2020-07-24  bugfix, was condensing multiple tabs ahead of quoted field

use File::Basename;


# This function was written to handle invalid CSV escaping in a consistent
# manner that will, hopefully, result in output fields that are closer to the
# intent of the malformed fields (guessing intent of malformed fields is,
# admittedly, fraught with peril).
#
# Due to the way that Excel handles quotes and using them to embed actual
# newlines (one line wrapped onto multiple lines), as well as how this script
# handles invalid CSV fields, the resulting tab-delimited text file may not
# load into Excel the same as if the .csv file had been loaded directly into
# Excel.  Thus the _not_excel part of the function name.  The resulting tsv
# file could, depending on how its quotes are arranged, result in strange
# Excel behavior.  The output should be fine for most other programs that read
# tab-delimited text, which is mainly what I wrote this for.  If you need
# something with better Excel compatability, I suggest using Python's
# csv.reader and csv.writer set to "dialect=csv.excel_tab".
#
#
# The following requirements of RFC 4180 are intentionally ignored:
#
# 1) "Fields containing line breaks (CRLF), double quotes, and commas
#     should be enclosed in double-quotes."
# 2) "Each record is located on a separate line, delimited by a line
#     break (CRLF)."
# 3) "Each line should contain the same number of fields throughout the file."
# 4) "The last field in the record must not be followed by a comma."
#
# I'm going to assume that all four of these are likely defined in the RFC
# due to the allowance of embedded newlines in 1).  The common implementation
# of 1) appears to assumes MS-DOS CRLF end-of-line (EOL) characters, rather
# than unix LF or a mix-and-match of some lines ending in CRLF and others in
# LF (frequently occurs after editing a CRLF file in Unix).  For example,
# Excel appears to assume that LF will be used to denote an embedded newline,
# with CRLF indicating a "real" newline.  This, as well as 2), is problematic
# in Unix software, which generally assumes LF is the EOL character, since
# lines containing embedded newlines are then split across multiple lines of
# text.  This is where 3) and 4) come in.  If we know that all lines contain
# the same number of fields, whenever we encounter a line that is too short,
# and its last field begins with a " that doesn't have a matching
# right-enclosing ", we can assume it is line-wrapped, read in the next line,
# and append it to the end of the original line, continuing to do so until we
# reach our fixed total number of fields per line.  Due to splitting of
# embedded newlines into multiple lines, escaped commas within embedded quotes
# could be left dangling at the end of a wrapped line instead of indicating a
# new field, so 4) is probably there to better handle these resulting
# problematic line-wrapped lines.  By requiring 3) and 4), we can sanity check
# the un-wrapped lines as well.  I wouldn't be surprised if 3) and 4) also
# make the embedded newline-aware parsing simpler?
#
# This tool is intended to be used within the Unix paradigm of line-based text
# processing tools, in which each line is generally independent of all other
# lines.  Allowing embedded newlines in 1) would break grep, sed, sort, wc,
# etc., and probably most Unix-based tools in general.  Wrapped lines could be
# reordered or filtered so that they are no longer sequential with their
# original line-wrapped partner lines or ordered for correct unwrapping.
#
# This function now supports embedded newlines as a non-default option, and
# should handle them similarly to Excel/Python on correctly quoted/escaped
# fields.  However, the output will strip the wrapped enclosing quotes as best
# it can, which may cause unexpected behavior in Excel if any malformed fields
# (lone double quotes) are left.  The output is meant to be read by software
# that doesn't support embedded newlines, but should still be parsed properly
# by Excel/Python as long as there aren't any malformed fields (lone double
# quotes) left to confuse them.  There is now also an option to escape double
# quotes and line-wrapped fields for improved Excel/Python compatability.
# Default settings replace embedded EOL with blanks or whitespace as
# appropriate to prevent non- whitespace from abutting, but can be set to
# output escaped (\\\r, \\\n) EOL instead.
#
# CRLF in 2) is not enforced, since we expect to receive data from Unix text
# files, which end in LF instead of CRLF.  Since 2), 3), and 4) are no longer
# needed to support 1), we have no reason to enforce them.  It is not uncommon
# for real data to have lines with variable numbers of fields.  After
# disallowing embedded newlines and ignoring 1, 2), and 3), we see no reason
# why terminal empty fields should be encoded as ,"" to satisfy 4) instead of
# simply using a terminal comma.  Excel also violates 4).  Thus, 4) is ignored
# as well, since it is frequently violated in real data.
#
# Additionally, the following requirement of RFC 4180 is loosened so as to
# attempt to accept invalid syntax wherever valid syntax could not be parsed.
# Excel also selectively loosens this requirement when reading in files:
#   
#   "If fields are not enclosed with double quotes, then double quotes may not
#    appear inside the fields."
#
#
# The specific invalid edge cases I tested for involve quotes inside fields
# that *aren't* enclosed in double-quotes, as well as leading/trailing spaces
# before or after enclosing double-quotes.
#
#
#   RFC 4180 on double-quotes:
#
#     "If double-quotes are used to enclose fields, then a double-quote
#      appearing inside a field must be escaped by preceding it with
#      another double quote.
#
#     "If fields are not enclosed with double quotes, then double quotes
#      may not appear inside the fields."
#
#
# This would appear to rule out treating fields surrounded in double-quotes
# as "enclosed" if there are any leading/trailing spaces, or if there are
# an even number of leading or trailing double-quotes.  However, this all
# depends on the defintion of "enclosed".
#
#
#   From Wikipedia's discussion on leading/trailing spaces:
#
#     "According to RFC 4180, spaces outside quotes in a field are not
#      allowed; however, the RFC also says that "Spaces are considered
#      part of a field and should not be ignored." and "Implementors
#      should 'be conservative in what you do, be liberal in what you
#      accept from others' (RFC 793, section 2.10) when processing CSV
#      files.""
#
#
# I choose to ignore leading/trailing spaces for the purposes of determining
# whether a field is enclosed within double-quotes or not.  I shall peserve
# them if the field is not enclosed ('spaces are considered part of a field
# and should not be ignored' and 'be conservative in what you do'), but remove
# them if the field is considered enclosed ('be liberal in what you accept
# from others').  Detecting and removing spaces that surround enclosing
# double-quotes allows for using spaces to visually align columns of variable
# width fields when displayed on a simple terminal or text editor that would,
# otherwise, not be vertically aligned.  I believe that this strikes a useful
# balance between simply preserving or removing all leading/trailing spaces.
# This attempts to follow the principle of "least surprise" to the user.
#
#
# This implementation results in the following behavior:
#
#   1) One or more embedded tab characters (not escaped as \t) in a row are
#      either removed (at beginning/end of field or space on either side),
#      or replaced with a single space (if removing would otherwise merge
#      non-space text together).  Existing spaces are preserved, only
#      one or more unescaped tabs in a row are stripped/condensed.
#
#      While not exactly standards compliant (is there even a standard for
#      for how to handle embedded tabs??), I have observed this to be closest
#      to original intent in the wild, where embedded tabs are almost always
#      due to user input error (such as pressing the tab key in an attempt to
#      advance to the next field/line in a form, etc.) and improper validation
#      or sanitizing of the data prior to storing it in a database.
#
#      Internal \r and \n that are not part of the EOL are treated the same
#      way as we treat embedded tabs, since they usually arise from similar
#      issues.  This behavior can be disabled by commenting out a single
#      substitution line in the function below.
#
#      Embedded newlines are, optionally, unwrapped and escaped, removing
#      their enclosed quotes similarly to as if Excel had parsed them.
#      This may trigger downstream false-positive line wraps in Excel.
#
#      Feel free to modify the handling of embedded tabs at the bottom of the
#      function below if you desire different tab handling behavior.
#
#   2) Leading/trailing spaces are ignored for the purposes of determining
#      if a field is enclosed in double-quotes.
#
#   3) "" is treated as escaped double-quotes, regardless of the position
#      within the field.  This includes the beginning/end of the field (with
#      or without leading/trailing spaces).  Therefore, in order to be
#      considered enclosed, the beginning/end of the field (after ignoring
#      leading/trailing spaces) must have an odd number of double-quotes,
#      so that the outer double-quote on each end is not treated as escaped.
#
#        !! ""aa,bb"" is not considered to be enclosed in double-quotes !!
#
#   4) Fields consisting solely of "", with or without leading/trailing
#      spaces, are special-cased as empty text.  I have now observed this
#      empty field double-quoting behavior in the wild, where it was clearly
#      meant to indicate an empty field.
#
#   5) Enclosing double-quotes are stripped from the output fields, along
#      with any leading/trailing spaces surrounding them.  Unenclosed
#      leading/trailing spaces are left as-is.
#
#   6) After removing enclosing double-quotes, "" are unescaped back into ",
#      regardless of where they are located within the field.
#
#   7) *MOST*, but not all, white space is preserved.  See 1) and 5).
#
#
# A relatively simple example: ""word""
#
# If we consider it to be enclosed within double-quotes, then the remaining
# single double-quotes are invalid, but if it is not considered enclosed, then
# the original "" are also invalid.  How best should we deal with it?
#
# Various tools result in different results:
#
#   word""    python and Excel, this looks not-right to me
#   "wor      many conversion websites, truncation is just plain wrong
#   "word"    this function, which I think handles it the most reasonably
#

{
  # wrap the sub in outer {} to implement static variables local to the {}
  # "state" declarations require perl >= v5.10.0,
  #  so the outer {} hack is more compatible
  $open_quotes_flag = 0;
  $prev_eol  = '';    # EOL at end of previous line wrap line
  $prev_estr = '';    # x or space depending on end of growing line unwrap

  sub csv2tsv_not_excel
  {
    my $line    = $_[0];
    my $opt_lw  = $_[1];    # detect and unwrap input line wraps
    my $opt_edq = $_[2];    # escape output double quotes
    my $delim   = $_[3];    # delimiter to use instead of comma
    my $opt_eol = $_[4];    # escape EOL instead of strip/whitespace
    my $c;
    my $i;
    my $j;
    my $k;
    my $n;
    my $strip_flag;
    my $eol;
    my @temp_array;        # only used for splitting on unenclosed fields
    my @temp_array2;
    my $quote_count;
    my $prepend_line;
    my $i_first_quotes;    # -1 if not line-wrapped
    my $i_last_quotes;     # @temp_array + 1 if not line-wrapped
    my $len1;
    my $len2;
    my $i_start;
    my $i_end;

    if (!defined($delim))
    {
        $delim = ',';
    }

    # placeholder strings unlikely to ever be encountered normally
    #    my $tab = '___TaBChaR___';
    #    my $dq  = '___DqUOteS___';
    #
    my $tab = "\x19";    # need single char for split regex
    my $dq  = "\x16";    # need single char for split regex
    
    # remove null characters, since I've only seen them randomly introduced
    #  during NTFS glitches; "Null characters may be inserted into or removed
    #  from a stream of data without affecting the information content of that
    #  stream."
    $line =~ s/\x00//g;
    
    # remove UTF8 byte order mark, since it corrupts the first field
    # also remove some weird corruption of the UTF8 byte order mark (??)
    #
    # $line =~ s/^\xEF\xBB\xBF//;      # actual BOM
    # $line =~ s/^\xEF\x3E\x3E\xBF//;  # corrupted BOM I have seen in the wild
    $line =~ s/^(\xEF\xBB\xBF|\xEF\x3E\x3E\xBF)//;

    # replace any (incredibly unlikely) instances of $dq with $tab
    $line =~ s/$dq/$tab/g;
    
    # replace embedded tabs with placeholder string, to deal with better later
    $line =~ s/\t/$tab/g;
    
    # HACK -- handle internal \r and \n the same way we handle tabs
    if ($opt_eol == 0)
    {
        # \v is vertical spacing, includes \r\n and other unicode:
        #   [\n\cK\f\r\x85\x{2028}\x{2029}]
        $line =~ s/[\v]+(?!$)/$tab/g;
    }
    # escape embedded EOL
    else
    {
        $line =~ s/\r(?!$)/\\\\\\r/g;
        $line =~ s/\n(?!$)/\\\\\\n/g;

        # replace remaining EOL-like characters as per usual
        # \v is vertical spacing, includes \r\n and other unicode:
        #   [\n\cK\f\r\x85\x{2028}\x{2029}]
        $line =~ s/[\v]+(?!$)/$tab/g;
    }
    
    # further escape ""
    $line =~ s/\"\"/$dq/g;

    # line-wrapped, delayed insert of escaped EOL from previous line
    # delay in order to not add escaped EOL to final line
    $prepend_line = '';
    if ($open_quotes_flag)
    {
        # replace EOL with space if growing text would abut
        # leave $dq out of the checks intentionally
        if ($opt_eol == 0 &&
            $prev_estr =~ /\S$/ && $line =~ /^\S/ &&
            !($line =~ /^[ $tab]*\"[ $tab]*(\Q$delim\E|[\r\n]*$)/))
        {
            $prepend_line .= " ";
            
            # line end book keeping string now ends in a space
            $prev_estr = ' ';
        }

        if ($opt_eol)
        {
            $prepend_line .= $prev_eol;
        }
    }

    # only apply slow tab expansion to lines still containing quotes
    @temp_array     = ();
    $i_last_quotes  = @temp_array + 1;
    $i_first_quotes = -1;
    
    if ($line =~ /\"/)
    {
        # slightly faster than split loop on rows with many quoted fields,
        #  but *much* slower on lines containing very few quoted fields
        # also cannot deal with line-wraps
        # use split loop instead
        #
        # /e to evaluate code to handle different capture cases correctly
        #$line =~ s/(\Q$delim\E?)((?<![^\Q$delim\E $tab$dq])"[^\t"]+"(?![^\Q$delim\E $tab$dq\r\n]))|(\Q$delim\E)/defined($3) ? "\t" : ((defined($1) && $1 ne '') ? "\t$2" : $2)/ge;

        # convert commas only if they are not within double quotes
        # incrementing $i within array access for minor speed increase
        #   requires initializing things to -2
        #@temp_array    = split /((?<![^\Q$delim\E $tab$dq])"[^\t"]+"(?![^\Q$delim\E $tab$dq\r\n]))/, $line, -1;
        #@temp_array    = split /((?:(?<=^)|(?<=\Q$delim\E))[ $tab$dq]*\"[^\t"]+\"[ $tab$dq\r\n]*(?=(?:\Q$delim\E|$)))/, $line, -1;
        #@temp_array    = split /((?:^|(?<=\Q$delim\E))[ $tab$dq]*\"[^\t"]+\"[ $tab$dq\r\n]*(?:$|(?=\Q$delim\E)))/, $line, -1;
        @temp_array    = split /((?:^|(?<=\Q$delim\E))[ $tab]*[$dq]*\"[^\t"]+\"[$dq]*[ $tab]*[\r\n]*(?:$|(?=\Q$delim\E)))/, $line, -1;
        $n             = @temp_array - 2;
        $i_last_quotes = @temp_array + 1;
        
        if ($opt_lw)
        {
          # scan for first quotes
          if ($open_quotes_flag)
          {
            for ($i = -1; $i < @temp_array - 1;)
            {
                if ($temp_array[$i += 1] =~ /\"/)
                {
                    $i_first_quotes = $i;

                    # these quotes close an existing line wrap
                    # reset the line end book keeping string
                    $prev_estr = '';

                    last;
                }
            }
          }

          # scan for potential line wrap
          for ($i = @temp_array - 1; $i >= 0; $i -= 2)
          {
            # abort, field is at or before closing quotes from earlier line
            if ($open_quotes_flag && $i < $i_first_quotes)
            {
                last;
            }
          
            # at or before field enclosed in quotes, abort line wrap check
            if ($i < @temp_array - 1 &&
                $temp_array[$i+1] =~ /\"/)
            {
                last;
            }
            
            if ($temp_array[$i] =~ /\"/)
            {
                # only accept valid escaped end fields
                if ($temp_array[$i] =~ /(^|(?:\Q$delim\E))[ $tab]*[$dq]*\"[^\"]*$/)
                {
                    # only a single " on entire line
                    # will be used to close a line-wrap; end check
                    @temp_array2 = split /\"/, $temp_array[$i], -1;
                    if ($i <= $i_first_quotes && @temp_array2 <= 2)
                    {
                        # either " closes the field and nothing comes
                        #  afterwards (@temp_array2 <= 1)
                        # or " followed by EOL character (@temp_array2 == 2)
                        if (@temp_array2 <= 1 ||
                            (@temp_array2 == 2 &&
                             $temp_array2[1] =~ /^[$dq]*[ $tab]*[\r\n]*$/))
                        {
                            last;
                        }
                    }

                    $i_last_quotes = $i;
                }

                # stop scaning after first encountered bare "
                last;
            }
          }
        }

        # convert all commas to tabs, no special line-wrap cases
        # this avoids costly checks inside the loop
        if ($open_quotes_flag == 0 && $i_last_quotes > @temp_array)
        {
            if ($delim eq ',')
            {
                for ($i = -2; $i < $n;)
                {
                    $temp_array[$i += 2] =~ tr/\,/\t/;
                }
            }
            elsif ($delim eq ';')
            {
                for ($i = -2; $i < $n;)
                {
                    $temp_array[$i += 2] =~ tr/\;/\t/;
                }
            }
            else
            {
                for ($i = -2; $i < $n;)
                {
                    $temp_array[$i += 2] =~ s/\Q$delim\E/\t/g;
                }
            }

            # reset the line end book keeping string
            # should have no impact on output, reset out of paranoia
            $prev_estr = '';
        }
        # deal with line wraps
        # lines within line-wraps, but without ", passed through as-is
        elsif ($line =~ /\"/)
        {
            $strip_flag = 1;

            for ($i = 0; $i < @temp_array; $i += 2)
            {
                # check for valid quotes in adjacent field
                if ($i < $i_first_quotes && $temp_array[$i+1] =~ /\"/)
                {
                    @temp_array2 = split /\"/, $temp_array[$i+1], -1;

                    # stop protecting commas after first " field,
                    # whether it was validly escaped or not
                    $open_quotes_flag = 0;
                    $strip_flag       = 0;

                    # convert the commas after the first quote
                    if ($delim eq ',')
                    {
                        for ($j = 1; $j < @temp_array2; $j++)
                        {
                            $temp_array2[$j] =~ tr/\,/\t/;
                        }
                    }
                    elsif ($delim eq ';')
                    {
                        for ($j = 1; $j < @temp_array2; $j++)
                        {
                            $temp_array2[$j] =~ tr/\;/\t/;
                        }
                    }
                    else
                    {
                        for ($j = 1; $j < @temp_array2; $j++)
                        {
                            $temp_array2[$j] =~ s/\Q$delim\E/\t/g;
                        }
                    }

                    $temp_array[$i+1] = join '"', @temp_array2;
                    
                    # remove the first quote
                    $temp_array[$i+1] =~ s/\"//;
                    
                    next;
                }

                # convert any fields not within line-wrap boundaries
                if ($open_quotes_flag == 0 &&
                    $i > $i_first_quotes && $i < $i_last_quotes)
                {
                    if ($delim eq ',')
                    {
                        $temp_array[$i] =~ tr/\,/\t/;
                    }
                    elsif ($delim eq ';')
                    {
                        $temp_array[$i] =~ tr/\;/\t/;
                    }
                    else
                    {
                        $temp_array[$i] =~ s/\Q$delim\E/\t/g;
                    }
                }
                # convert commas within line-wrap boundaries
                elsif ($temp_array[$i] =~ /\"/)
                {
                    @temp_array2 = split /\"/, $temp_array[$i], -1;

                    # stop protecting commas after first " field,
                    # whether it was validly escaped or not
                    if ($open_quotes_flag)
                    {
                        $open_quotes_flag = 0;
                        $strip_flag       = 1;
                    }
                    else
                    {
                        # begin converting up to new opening quotes
                        if ($i <= $i_last_quotes)
                        {
                            if ($i < $i_last_quotes || 0 < @temp_array2 - 1)
                            {
                                if ($delim eq ',')
                                {
                                    $temp_array2[0] =~ tr/\,/\t/;
                                }
                                elsif ($delim eq ';')
                                {
                                    $temp_array2[0] =~ tr/\;/\t/;
                                }
                                else
                                {
                                    $temp_array2[0] =~ s/\Q$delim\E/\t/g;
                                }
                            }
                        }
                    }

                    # continue converting, handle new opening quotes
                    for ($j = 1; $j < @temp_array2; $j++)
                    {
                        # open new line wrap
                        if ($i == $i_last_quotes && $j == @temp_array2 - 1)
                        {
                            # reset the line end book keeping string
                            # should already be reset, but reset just in case
                            $prev_estr = '';

                            # remove opening " from line-wrap
                            if (@temp_array2 > 1)
                            {
                                $k = @temp_array2 - 2;
                                
                                # remove padding spaces before opening "
                                $temp_array2[$k] =~ s/[ $tab]+([$dq]*)$/$1/;
                                
                                # remove opening " by tacking on remaining,
                                # then deleting remaining
                                $temp_array2[$k] .= $temp_array2[$k+1];
                                delete $temp_array2[$k+1];
                            }
                            $temp_array[$i] = join '"', @temp_array2;
                        }
                        # convert commas to tabs
                        elsif ($i < $i_last_quotes || $j < @temp_array2 - 1)
                        {
                            if ($delim eq ',')
                            {
                                $temp_array2[$j] =~ tr/\,/\t/;
                            }
                            elsif ($delim eq ';')
                            {
                                $temp_array2[$j] =~ tr/\;/\t/;
                            }
                            if ($delim eq ',')
                            {
                                $temp_array2[$j] =~ s/\Q$delim\E/\t/g;
                            }
                        }
                    }

                    $temp_array[$i] = join '"', @temp_array2;
                }
            }
            
            # remove proper line wrap terminating enclosing quotes
            if ($strip_flag)
            {
                $temp_array[0] =~ s/\"([$dq]*)[ $tab]*((?:\t|[\r\n]*$))/$1$2/;
            }
        }

        # set next line to deal with line-wrap
        if ($i_last_quotes <= @temp_array)
        {
            $open_quotes_flag = 1;
        }
        
        $line = join '', @temp_array;
    }
    elsif ($open_quotes_flag == 0)
    {
        if ($delim eq ',')
        {
            $line =~ tr/\,/\t/;
        }
        elsif ($delim eq ';')
        {
            $line =~ tr/\;/\t/;
        }
        else
        {
            $line =~ s/\Q$delim\E/\t/g;
        }

        # reset the line end book keeping string
        # should have no impact on output, reset out of paranoia
        $prev_estr = '';
    }
    
    # escape EOL if line-wrapped
    if ($open_quotes_flag)
    {
        # strip EOL, escape it and save it for potential future use later
        $prev_eol = '';
        if ($line =~ s/([\r\n]+)$//)
        {
            $prev_eol = $1;
            $prev_eol =~ s/\r/\\\\\\r/g;
            $prev_eol =~ s/\n/\\\\\\n/g;
        }
        
        # line end book keeping string now ends in a space
        if ($line =~ /[ $tab]$/)
        {
            $prev_estr = " ";
        }

        # line end book keeping string now ends in not-whitespace
        if ($line =~ /\S$/)
        {
            $prev_estr = "x";
        }
        
        # lines ending in control character whitespace don't alter book keeping
    }
    
    # finish dealing with embedded tabs
    # remove tabs entirely, preserving surrounding whitespace
    $line =~ s/(\s|^)($tab)+/$1/g;
    $line =~ s/($tab)+(\s|$)/$2/g;
    # replace remaining tabs with spaces so text doesn't abut together
    $line =~ s/($tab)+/ /g;

    # clean up spaces around quotes, handle single "" on entire line
    # only apply to non- line wrapped fields
    if ($line =~ /[\"$dq]/)
    {
      if ($open_quotes_flag ||
          $i_first_quotes >= 0 ||
          $i_last_quotes  <= @temp_array)
      {
        # unescape "" before space/quotes removal
        $line =~ s/$dq/\"\"/g;

        @temp_array2 = split /\t/, $line, -1;

        $eol = '';
        if ($temp_array2[@temp_array2 - 1] =~ s/([\r\n]+)$//)
        {
            $eol = $1;
        }
        

        # skip when entire line is part of a line wrap
        $i_start = 0;
        $i_end   = @temp_array2;
        if ($open_quotes_flag && @temp_array2 == 1)
        {
            $i_end = -1;
        }
        # skip only the portion within a line wrap
        else
        {
            # skip first field of newly ended line wrap
            if ($i_first_quotes >= 0)
            {
                $i_start = 1;
            }

            # skip last field of newly opened line wrap
            if ($i_last_quotes <= @temp_array)
            {
                $i_end--;
            }
        }
        

        for ($i = $i_start; $i < $i_end; $i++)
        {
            # find start run of quotes
            if ($temp_array2[$i] =~ /^ *(\"+)/)
            {
                $len1 = length $1;
                
                # find end run of quotes
                #
                # slightly faster to reverse then search forwards
                # than to search for match from end of string

                if (reverse($temp_array2[$i]) =~ /^ *(\"+)/)
                {
                    $len2 = length $1;

                    # remove enclosing spaces and quotes
                    # considered enclosed if both start/end are odd
                    if (($len1 & 1) && ($len2 & 1))
                    {
                        # .*? isn't needed here
                        # since we have an end-matching constraint
                        $temp_array2[$i] =~ s/^ *\"(.*)\" *$/$1/;
                    }
                    # or if entire field is even number of quotes
                    # must be even at this point due to first check
                    else
                    {
                        # .*? isn't needed here
                        # since we have an end-matching constraint
                        $temp_array2[$i] =~ s/^ *\"(\"*)\" *$/$1/;
                    }
                }
            }
        }

        # add EOL back to end of line
        if ($eol ne '')
        {
            $temp_array2[@temp_array2 - 1] .= $eol;
        }
        
        $line = join "\t", @temp_array2;
      }

      # safe to strip entire line at once
      else
      {
        # remove enclosing spaces and quotes
        $line =~ s/(?:(?<=\t)|^) *([$dq]*)\"([^\t]*?)\" *(?=\t|[\r\n]*$)/$1$2/g;

        # unescape "" after first pass of space/quotes removal
        $line =~ s/$dq/\"\"/g;

        # remove enclosing spaces and quotes from all-quotes fields
        $line =~ s/(?:(?<=\t)|^) *\"(\"*?)\" *(?=\t|[\r\n]*$)/$1/g;
      }
    }

    # unescape escaped double-quotes
    $line =~ s/\"\"/\"/g;

    # Escape all final quotes to prevent Excel from misinterpreting them.
    #
    # Escape first/last fields of line wraps, since figuring out if there
    # are any final double quotes within multiple lines of wrapping isn't
    # simple.
    if ($opt_edq &&
        ($line           =~ /\"/ ||
         $i_first_quotes >= 0    ||
         $i_last_quotes  <= @temp_array))
    {
        # remove the terminal EOL character(s)
        $eol = '';
        if ($line =~ s/([\r\n]+)$//)
        {
            $eol = $1;
        }
        
        @temp_array2 = split /\t/, $line, -1;
        for ($i = 0; $i < @temp_array2; $i++)
        {
            # first field of newly ended line wrap
            if ($i_first_quotes >= 0 && $i == 0)
            {
                $temp_array2[$i] =~ s/\"/\"\"/g;
                $temp_array2[$i] = $temp_array2[$i] . '"';
            }
            # last field of newly opened line wrap
            elsif ($i_last_quotes <= @temp_array && $i == @temp_array2 - 1)
            {
                $temp_array2[$i] =~ s/\"/\"\"/g;
                $temp_array2[$i] = '"' . $temp_array2[$i];
            }
            # entire line is within a line wrap
            elsif ($open_quotes_flag && @temp_array2 == 1)
            {
                # only escape the double quotes, don't enclose them
                $temp_array2[$i] =~ s/\"/\"\"/g;
            }
            # all other fields with quotes in them
            elsif ($temp_array2[$i] =~ s/\"/\"\"/g)
            {
                $temp_array2[$i] = '"' . $temp_array2[$i] . '"';
            }
        }
        
        # closing lone " on a line by itself
        if (@temp_array2 == 0 && $i_first_quotes >= 0)
        {
            $temp_array2[0] = '"';
        }
        # opening lone " on a line by itself
        if (@temp_array2 == 0 && $i_last_quotes <= @temp_array)
        {
            $temp_array2[0] = '"';
        }
        
        # instantiate $temp_array2[0] if it doesn't exist,
        # so that adding the EOL won't write out of bounds
        if (@temp_array2 == 0)
        {
            $temp_array2[0] = '';
        }
        
        # add the EOL back in after enclosing in ""
        #
        # should be faster to add EOL to last field before joining
        # instead of adding it to the end of the joined line
        $temp_array2[@temp_array2 - 1] .= $eol;

        $line = join "\t", @temp_array2;
    }
    
    # prepend any line wrap related characters from previous line(s)
    $line = $prepend_line . $line;
    
    return $line;
  }
}



# begin main()

$opt_delim      = ',';
$opt_line_wrap  = 0;
$opt_escape_dq  = 0;
$opt_escape_eol = 0;

$error_flag     = 0;

# read in command line arguments
$num_files = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--unwrap/i)
        {
            $opt_line_wrap = 1;
        }
        elsif ($field =~ /^--escape-dq/i)
        {
            $opt_escape_dq = 1;
        }
        elsif ($field =~ /^--escape-eol/i)
        {
            $opt_escape_eol = 1;
        }
        elsif ($field =~ /^--semicolon/i)
        {
            $opt_delim = ';';
        }
        else
        {
            $error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $filename = $field;
            $num_files++;
        }
    }
}


if ($error_flag)
{
    $program_name = basename($0);
    
    printf STDERR "Usage: $program_name [options] input.csv\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --escape-dq    escape \" for better Excel compatability\n";
    printf STDERR "    --escape-eol   escape embedded EOL with \\\\\\\n";
    printf STDERR "    --semicolon    use semicolon as input delimiter instead of comma\n";
    printf STDERR "    --unwrap       unwrap lines with embedded newlines\n";
    
    exit(1);
}


# default to stdin if no filename given
if ($num_files == 0)
{
    $filename = '-';
    $num_files = 1;
}


open FILE, "$filename" or die "ABORT -- can't open file $filename\n";

while(defined($line=<FILE>))
{
    $line_new = csv2tsv_not_excel($line, $opt_line_wrap, $opt_escape_dq,
                                  $opt_delim, $opt_escape_eol);
    
    print "$line_new";
    
}
close FILE;

# print terminal EOL if the line-wrap was never closed
if ($open_quotes_flag)
{
    print "\n";
}
