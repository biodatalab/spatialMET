#!/usr/bin/perl -w

# make sure we treat the input file as ut8f
# this also screws things up for some reason, so don't do it
#use open ':encoding(utf8)';

# Don't set the output to UTF8.
# For some reason, it corrupt UTF8 characters into who knows what on output.
# If we set binmode utf8, we need to call utf8::decode() before printing
#  each line.  I still don't understand why things work this way.
# 
#binmode(STDOUT,':utf8');

# | can occur in metabolite names (!)
# so, use || as the delimiter
$delim_out = '||';

# regular characters, de-smartify them
$html_char_map_hash{'&amp;'}    = '&';
$html_char_map_hash{'&gt;'}     = '>';
$html_char_map_hash{'&lt;'}     = '<';
$html_char_map_hash{'&rsquo;'}  = "'";
$html_char_map_hash{'&lsquo;'}  = "'";
$html_char_map_hash{'&prime;'}  = "'";
$html_char_map_hash{'&ndash;'}  = '-';          # dash in number ranges

# UTF8
if (1)
{
    $html_char_map_hash{'&Delta;'}  = "\xCE\x94";
    $html_char_map_hash{'&alpha;'}  = "\xCE\xB1";
    $html_char_map_hash{'&beta;'}   = "\xCE\xB2";
    $html_char_map_hash{'&gamma;'}  = "\xCE\xB3";
    $html_char_map_hash{'&delta;'}  = "\xCE\xB4";
    $html_char_map_hash{'&mu;'}     = "\xCE\xBC";
    $html_char_map_hash{'&omega;'}  = "\xCF\x89";

    $html_char_map_hash{'&reg;'}    = "\xC2\xAE";	# registered (R)
    $html_char_map_hash{'&deg;'}    = "\xC2\xB0";	# degree
    $html_char_map_hash{'&szlig;'}  = "\xC3\x9F";	# sharp S, German B-like symbol
    $html_char_map_hash{'&ntilde;'} = "\xC3\xB1";	# n with tilde over it

    $html_char_map_hash{'&bull;'}   = "\xE2\x80\xA2";	# bullet
    $html_char_map_hash{'&trade;'}  = "\xE2\x84\xA2";	# trademark (TM)
}
# convert to plain ASCII text
else
{
    $html_char_map_hash{'&Delta;'}  = "delta ";
    $html_char_map_hash{'&alpha;'}  = "alpha";
    $html_char_map_hash{'&beta;'}   = "beta";
    $html_char_map_hash{'&bull;'}   = "* ";	# bullet
    $html_char_map_hash{'&deg;'}    = " degrees ";	# degree
    $html_char_map_hash{'&delta;'}  = "delta";
    $html_char_map_hash{'&gamma;'}  = "gamma";
    $html_char_map_hash{'&mu;'}     = "mu";
    $html_char_map_hash{'&ntilde;'} = "n";	# n with tilde over it
    $html_char_map_hash{'&omega;'}  = "omega";
    $html_char_map_hash{'&reg;'}    = "(R)";	# registered (R)
    $html_char_map_hash{'&szlig;'}  = "sz";	# sharp S, German B-like symbol
    $html_char_map_hash{'&trade;'}  = "(TM)";	# trademark (TM)
}

# debug HTML --> Windows ANSI characters
if (0)
{
    foreach $html (sort keys %html_char_map_hash)
    {
        $char = $html_char_map_hash{$html};
        #utf8::decode($char);
    
        if ($char ne '')
        {
            printf "%s\t%s\n", $html, $char;
        }
    }
    die;
}


$filename = shift;

open INFILE, "$filename" or die "ABORT -- can't open $filename\n";

$row         = -1;   # which 0-level indent we're currently on
$level       = -1;   # level of indention
@level_array = ();   # individual header at each level

%seen_long_header_hash = ();

$line_num = -1;
while(defined($line=<INFILE>))
{
    $line_num++;

    $line      =~ s/[\r\n]+//g;
    
    # skip initial XML header line
    # skip hmdb lines -- we know all the fields come from HMDB already
    if ($line =~ /^<(\?xml|hmdb xmlns|\/hmdb>)/)
    {
        next;
    }
    
    # remove leading/trailing spaces from line
    $line =~ s/(^\s+|\s+$)//;

    # back up line for debug printing
    $line_old = $line;


  # parse line until we have no non-whitespace left on it
  while ($line =~ /\S/)
  {
    # start of new field
    if ($line =~ m/^<(?!\/)([^<]+)>/)
    {
        $field        = $1;

        # remove parsed portion from beginning of line
        $rest_of_line = substr $line, (length $field) + 2;
        $line         = $rest_of_line;
        
        # evidently, if the end of the field is /
        # then it means the field does not exist, and there is no closing </>
        
        # expect a closing </> later
        if (!($field =~ /\/$/))
        {
            # start of new row
            if ($level == -1)
            {
                $row++;

                printf STDERR "Row:\t%d\n", $row;
            }

            # increment level
            $level++;
            $level_array[$level] = $1;
        }
        # else do nothing, due to ending in /
        
        $data = '';
    }

    $c = substr $line, 0, 1;
    
    # data to store for a field
    if ($line =~ /\S/ && $c ne '<')
    {
        # temporarily strip next closing </> and beyond from line
        $temp = $line;
        $temp =~ s/<\/[^<]+>.*//;
        
        # check for potential whitespace issues between lines
        if ($data ne '' &&
            !($data =~ /\s$/))
        {
            printf STDERR "WARNING: insert whitespace between wrapped lines: %d\n",
                $line_num;
        
            $data .= ' ';
        }
        
        $data .= $temp;

        # remove parsed portion from beginning of line
        $rest_of_line = substr $line, (length $temp);
        $line         = $rest_of_line;
    }
    
    # close a field
    if ($line =~ m/^<\/([^<]+)>/)
    {
        $skip_flag = 0;
        $field     = $1;

        if ($data =~ /\S/)
        {
            # strip tabs from end of field, since they are most
            # commonly data entry errors (someone hit tab to advance to
            # the next field and it was captured as text).
            $data =~ s/\t+$//;
            
            # strip any remaining whitespace at the end of the field
            $data =~ s/\s+$//;
            
            # there probably shouldn't be any whitespace at the beginning
            # of the field either, so make sure we strip it too
            $data =~ s/^\s+//;

            # deal with any remaining embedded tabs
            # remove tabs entirely, preserving surrounding spaces
            $data =~ s/( |^)(\t)+/$1/g;
            $data =~ s/(\t)+( |$)/$2/g;
            # replace remaining tabs with spaces so text doesn't abutt together
            $data =~ s/(\t)+/ /g;
        
            # deal with HTML-ized characters
            while ($data =~ m/(\&[A-Za-z]+;)/g)
            {
                $html_char = $1;
                
                # map it to its regular character
                $char_mapped = $html_char_map_hash{$html_char};
                if (defined($char_mapped))
                {
                    $data =~ s/\Q$html_char\E/$char_mapped/g;
                }
                # not hard coded in our list, store it for later
                else
                {
                    $seen_html_char_hash{$html_char} = 1;
                }
            }
        }
        
        # HACK -- <term> is used for flexible ontology headers
        # replace <term> with the value of data
        $term_flag = 0;
        if ($level > 0 && $level_array[1] eq 'ontology' &&
            $field eq 'term')
        {
            $level_array[$level - 1] = $data;

            $term_flag = 1;
        }

        # HACK -- <kind> is used for flexible ontology headers
        # replace <kind> with the value of data
        $kind_flag = 0;
        if ($level > 0 &&
            $level_array[1] =~ /^(experimental|predicted)_properties$/ &&
            $field eq 'kind')
        {
            $level_array[$level - 1] = $data;

            $kind_flag = 1;
        }
        
        # HACK -- <biospecimen> is used for flexible ontology headers
        # replace <biospecimen with the value of data
        $abconc_flag = 0;
        if ($level > 0 &&
            $level_array[1] eq 'abnormal_concentrations' &&
            $field eq 'biospecimen')
        {
            $level_array[$level - 1] = $data;

            $abconc_flag = 1;
        }
        

        # remove parsed portion from beginning of line
        $rest_of_line = substr $line, (length $field) + 3;
        $line         = $rest_of_line;
        
        $header = join "::", @level_array;

        # HACK -- condense terms
        if ($term_flag)
        {
            $header =~ s/::term$//;
            $header =~ s/::$data//;
        }
        if ($kind_flag)
        {
            $header =~ s/::kind$//;
            $header =~ s/::$data//;
        }
        if ($abconc_flag)
        {
            #$header =~ s/::biospecimen$//;
            $header =~ s/::$data//;
        }

        # all fields start with metabolite, so remove it
        $header =~ s/(?:^|\t)metabolite:://g;
        
        # creation dates/version of the database
        # might be the same for all metabolites
        if ($header eq 'creation_date' ||
            $header eq 'update_date' ||
            $header eq 'version')
        {
            $skip_flag = 1;
        }
        
        # there's just too much abnormal conecentration stuff
        # keep only the initial merged biospecimen
        if ($level > 0 &&
            $level_array[1] eq 'abnormal_concentrations')
        {
            if ($header =~ /^abnormal_concentrations::/ &&
                $header ne 'abnormal_concentrations::biospecimen')
            {
                $skip_flag = 1;
            }
        }

        # HACK -- more _properties clutter
        if ($level > 0 &&
            $level_array[1] =~ /^(experimental|predicted)_properties$/)
        {
            if ($header =~ /^(experimental|predicted)_properties$/ ||
                $field eq 'source')
            {
                $skip_flag = 1;
            }
            
            $header =~ s/::value$//;
        }

        # HACK -- more ontology clutter
        if ($level > 0 &&
            $level_array[1] eq 'ontology')
        {
            # remove oddly organized hierarchy stuff
            $header = $header;
            if ($header =~ /^ontology/)
            {
                $header =~ s/::descendants//g;
            }

            if ($header eq 'ontology' ||
                $field eq 'level' ||
                $field eq 'parent_id' ||
                $field eq 'type' ||
                $field eq 'definition' ||
                $field eq 'synonym')
            {
                $skip_flag = 1;
            }
            
            # condense more stuff
            if ($skip_flag == 0)
            {
                # condense Fungi, Microbes, and Plants
                $header =~
                  s/(::Source::Biological::(?:Fungi|Microbe|Plant)).*/$1/;

                # condense individual disease-related pathways
                $header =~ s/(::Biochemical pathway)::.*/$1::Disease/;

                # condense individual health conditions
                $header =~ s/(::Health condition)::.*/$1::Details/;
            }
        }
        
        if ($data =~ /\S/ && $skip_flag == 0)
        {
            # printf STDERR "%d\t%s\t%s\n", $row, $header, $data;
        
            $data_array_hash[$row]{$header}{$data} = 1;
            $seen_long_header_hash{$header}        = 1;
        }
        
        # remove element and back up a level
        # delete is discourage on array elements, so use splice instead
        if ($level >= 0)
        {
            splice @level_array, $level, 1;
        }
        # uh oh, this shouldn't happen
        else
        {
            $level = -1;

            printf STDERR "WARNING -- parsing error at line %d: %s\n",
                $line_num, $line_old;
        }

        $level--;
        
        $data = '';
    }
  }    
}
$num_rows = $row + 1;


# print header line
@header_array = sort keys %seen_long_header_hash;
print "RowIndex";
for ($i = 0; $i < @header_array; $i++)
{
    $header = $header_array[$i];

    print "\t";
    print $header;
}
print "\n";



# print data
for ($row = 0; $row < $num_rows; $row++)
{
    print "$row";

    for ($i = 0; $i < @header_array; $i++)
    {
        $data_str    = '';
        $header = $header_array[$i];

        @temp_array = ();
        if (defined($data_array_hash[$row]{$header}))
        {
            @temp_array = sort keys %{$data_array_hash[$row]{$header}};
        
            $data_str = join $delim_out, @temp_array;
        }
        
        $len_data_str = length $data_str;
        
        # too long for Excel
        # Excel documentation lies, the limit is empirically less than 32767
        #
        # NOTE -- only the first 1000 characters are displayed
        #
        $excel_width_limit = 32758;
        if ($len_data_str > $excel_width_limit)
        {
            $grow_str = '';
            $data_str = 'EXCEL_TRUNCATED';

            for ($j = 0; $j < @temp_array; $j++)
            {
                if ($j)
                {
                    $grow_str .= $delim_out . $temp_array[$j];
                }
                else
                {
                    $grow_str  = $temp_array[$j];
                }

                $test_str  = $grow_str . $delim_out . 'EXCEL_TRUNCATED';
                
                if (length $test_str <= $excel_width_limit)
                {
                    $data_str = $test_str;
                }
            }
            
            if (!defined($truncated_hash{$header}))
            {
                $truncated_hash{$header} = 0;
            }
            $truncated_hash{$header} += 1;
        }

        # scan data string for unicode characters
        $unicode = $data_str;
        $unicode =~ s/[\x00-\x7F]//g;
        
        if ($unicode ne '')
        {
            # first, decode the unicode string
            # into single characters, so substr works correctly
            utf8::decode($unicode);
            $len = length $unicode;
            
            for ($j = 0; $j < $len; $j++)
            {
                $c = substr $unicode, $j, 1;
            
                # encode it back again, since input is multi-byte chars
                utf8::encode($c);
                
                $hex_str = uc unpack("H*", $c);
                
                $seen_unicode_hash{$hex_str} = $c;
            }
        }
    
        print "\t";
        print $data_str;
    }
    print "\n";
}


# print list of headers that were truncated for Excel
@truncated_array = sort keys %truncated_hash;
foreach $header (@truncated_array)
{
    printf STDERR "Excel-truncated column:\t%d\t%s\n",
        $truncated_hash{$header}, $header;
}

# print list of encountered HTML-ized characters we still need to deal with
@html_array = sort keys %seen_html_char_hash;
foreach $html_char (@html_array)
{
    printf STDERR "HTML character:\t%s\n", $html_char;
}


# print list of encountered unicode characters
foreach $hex_str (sort keys %seen_unicode_hash)
{
    $c = $seen_unicode_hash{$hex_str};
    
    #printf STDERR "Unicode character:\t%s\t%s\n", $c, $hex_str;
}
