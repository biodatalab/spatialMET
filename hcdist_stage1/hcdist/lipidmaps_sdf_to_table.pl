#!/usr/bin/perl -w

# 2023-06-19:  print LM_ID in corrected semicolon spacing debug messages
# 2023-06-15:  detect and correct Cer(x; y/...) space after semicolon
# 2022-10-10:  initial commit to git


use File::Basename;


sub cmp_keys
{
    my $known_order_a = $output_known_header_order_hash{$a};
    my $known_order_b = $output_known_header_order_hash{$b};
    my $count_a = $key_count_hash{$a};
    my $count_b = $key_count_hash{$b};
    
    # most important annotation headers first
    if (defined($known_order_a) && !defined($known_order_b)) { return -1; }
    if (defined($known_order_b) && !defined($known_order_a)) { return  1; }

    # group most important headers manually
    if (defined($known_order_a) && defined($known_order_b))
    {
        return $known_order_a <=> $known_order_b;
    }
    
    # sort by number of values present
    if ($count_a > $count_b) { return -1; }
    if ($count_b > $count_a) { return  1; }
    
    # sort by consensus order within sdf file
    return ($key_rank_hash{$a} <=> $key_rank_hash{$b});
}


%seen_entry_row_hash = ();
%table_hash          = ();
$key_rank_hash       = ();    # sum rank of key, divide into avg later
$key_count_hash      = ();    # how many times we saw this key
$num_keys            =  0;

$i = 0;
$output_known_header_order_hash{'LM_ID'} = $i++;
$output_known_header_order_hash{'EXACT_MASS'} = $i++;
$output_known_header_order_hash{'FORMULA'} = $i++;
$output_known_header_order_hash{'CATEGORY'} = $i++;
$output_known_header_order_hash{'MAIN_CLASS'} = $i++;
$output_known_header_order_hash{'SUB_CLASS'} = $i++;
$output_known_header_order_hash{'CLASS_LEVEL4'} = $i++;
$output_known_header_order_hash{'CHOSEN_NAME'} = $i++;
$output_known_header_order_hash{'NAME'} = $i++;
$output_known_header_order_hash{'SYSTEMATIC_NAME'} = $i++;
$output_known_header_order_hash{'SYNONYMS'} = $i++;
$output_known_header_order_hash{'ABBREVIATION'} = $i++;
$output_known_header_order_hash{'ISO'} = $i++;
$output_known_header_order_hash{'INCHI'} = $i++;
$output_known_header_order_hash{'INCHI_KEY'} = $i++;
$output_known_header_order_hash{'SMILES'} = $i++;





$filename = shift;

# read from STDIN if no filename is specified
if (!defined($filename))
{
    $filename = '-';
}


if ($filename =~ /^-/ && length $filename >= 2)
{
    $program_name = basename($0);
    
    print "Usage: $program_name lipidmaps_sdf.sdf\n";
    exit(1);
}


open INFILE, "$filename" or die "ABORT -- cannot open file $filename\n";


$row       = 0;
$entry_row = 0;
$key_rank  = 0;
$lm_id     = '';
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    $c = substr $line, 0, 1;
    
    # first row of a new entry
    #
    # NOTE -- most entries start with a line containing the identifier
    #         but some have missing or duplicate identifiers,
    #         so have to use the row we started on as the identifier
    # 
    if ($key_rank == 0 && $c =~ /\S/)
    {
        $entry_row = $row;
        $seen_entry_row_hash{$entry_row} = 1;

        $key_rank = 1;
        $lm_id = '';
    }
    
    if ($c eq '>')
    {
        if ($line =~ /^> *<([^>]+)>/)
        {
            $key = $1;
            
            while(defined($line=<INFILE>) && $line =~ /\S/)
            {
                $line      =~ s/[\r\n]+//g;
                $line_orig = $line;
            
                # sigh, the usual crap happens in LIPID MAPS too...
                $line =~ s/^;*\s+;+//;
                $line =~ s/;*\s+;*$//;
                
                # convert / into ; to separate sub-fields
                if ($line =~ / \/ /)
                {
                    $line =~ s/\s+\/\s+/; /g;
                    
                    $key_with_slash_hash{$key} = 1;
                }
                
                
                # Correct for obvious spacing errors after ;
                #
                # Otherwise, it splits into two partial strings,
                # just like on the LipidMaps website (also errors).
                #
                if ($line =~ s/(\([^\(\)]+;)\s+([^\(\)]+:[^\(\)]+\/)/$1$2/)
                {
                    printf STDERR "SPACING_ERROR:\t%s\t%s\n",
                        $lm_id, $line_orig;
                    printf STDERR "SPACING_FIXED:\t%s\t%s\n",
                        $lm_id, $line;
                }
                
                
                # store current LM_ID for debug printing
                if ($key eq 'LM_ID')
                {
                    $lm_id = $line;
                }
                
                $table_hash{$entry_row}{$key}{$line} = 1;

                if (!defined($key_count_hash{$key}))
                {
                    $key_count_hash{$key} = 0;
                    $key_rank_hash{$key}  = 0;

                    $num_keys++;
                }

                $key_rank_hash{$key}  += $key_rank;
                $key_count_hash{$key} += 1;
            }
            
            $key_rank++;
        }
    }
    elsif ($c eq '$' && $line =~ /^\$\$\$\$/)
    {
        $key_rank = 0;
    }
    
    $row++;
}


# calculate average ranks for each key
foreach $key (%key_rank_hash)
{
    if ($key_count_hash{$key})
    {
        $key_rank_hash{$key} /= $key_count_hash{$key};
    }
}

@entry_row_array  = sort keys %table_hash;
@sorted_key_array = sort cmp_keys keys %key_rank_hash;


# strip [iso#] from fields, split into separate new field
foreach $entry_row (@entry_row_array)
{
    for ($i = 0; $i < @sorted_key_array; $i++)
    {
        $key = $sorted_key_array[$i];
        
        if (defined($table_hash{$entry_row}{$key}))
        {
            @value_str_array = sort keys %{$table_hash{$entry_row}{$key}};
            
            foreach $value_str (@value_str_array)
            {
                @match_array = $value_str =~ m/\[(iso[0-9]+)\]/g;
                
                $match_flag = 0;
                foreach $match (@match_array)
                {
                    $table_hash{$entry_row}{'ISO'}{$match} = 1;

                    $key_rank_hash{'ISO'}  += $key_rank;
                    $key_count_hash{'ISO'} += 1;
                    $match_flag = 1;
                }
                
                # replace value with value stripped of iso
                if ($match_flag)
                {
                    delete $table_hash{$entry_row}{$key}{$value_str};
                    $value_str =~ s/\s*\[(iso[0-9]+)\]\s*//g;
                    if ($value_str =~ /\S/)
                    {
                        $table_hash{$entry_row}{$key}{$value_str} = 1;
                    }
                }
            }
        }
    }
}

if (defined($key_rank_hash{'ISO'}))
{
    $key_rank_hash{'ISO'} /= $key_count_hash{'ISO'};
    @sorted_key_array = sort cmp_keys keys %key_rank_hash;
}

# LipidMap can have NAME present with blank SYSTEMATIC_NAME,
# as well as SYSTEMATIC_NAME present with no NAME.
#
# Generate a new column that prefers one over the other
if (defined($key_rank_hash{NAME}) &&
    defined($key_rank_hash{SYSTEMATIC_NAME}))
{
    $key = 'CHOSEN_NAME';

    if (!defined($key_count_hash{$key}))
    {
        $key_count_hash{$key} = 0;
        $key_rank_hash{$key}  = $num_keys;

        $num_keys++;
    }

    foreach $entry_row (@entry_row_array)
    {
        $name     = '';
        $sys_name = '';
        
        if (defined($table_hash{$entry_row}{NAME}))
        {
            @value_array = sort keys %{$table_hash{$entry_row}{NAME}};
            $value_str = join '; ', @value_array;
            
            $name = $value_str;
        }
        if (defined($table_hash{$entry_row}{SYSTEMATIC_NAME}))
        {
            @value_array = sort keys %{$table_hash{$entry_row}{SYSTEMATIC_NAME}};
            $value_str = join '; ', @value_array;
            
            $sys_name = $value_str;
        }

        $chosen_name = '';
        
        # pick shorter of the two names
        #
        # Example: LMST03020212
        #   NAME:             (10Z)-19-fluoro-1alpha,25-dihydroxyvitamin D3 / (10Z)-19-fluoro-1alpha,25-dihydroxycholecalciferol
        #   SYSTEMATIC_NAME:  (5Z,7E,10Z)-(1S,3R)-19-fluoro-9,10-seco-5,7,10(19)-cholestatriene-1,3,25-triol
        #        
        if ($name ne '' && $sys_name ne '')
        {
            $chosen_name = $name;
            if (length $sys_name < length $name)
            {
                $chosen_name = $sys_name;
            }
        }
        elsif ($name ne '')
        {
            $chosen_name = $name;
        }
        elsif ($sys_name ne '')
        {
            $chosen_name = $sys_name;
        }
        # get desperate
        else
        {
            $class_level4 = '';
            $sub_class    = '';
            $main_class   = '';
            $category     = '';
            
            if (defined($table_hash{$entry_row}{CLASS_LEVEL4}))
            {
                @value_array = sort keys %{$table_hash{$entry_row}{CLASS_LEVEL4}};
                $value_str = join '; ', @value_array;
            
                $class_level4 = $value_str;
            }
            if (defined($table_hash{$entry_row}{SUB_CLASS}))
            {
                @value_array = sort keys %{$table_hash{$entry_row}{SUB_CLASS}};
                $value_str = join '; ', @value_array;
            
                $sub_class = $value_str;
            }
            if (defined($table_hash{$entry_row}{MAIN_CLASS}))
            {
                @value_array = sort keys %{$table_hash{$entry_row}{MAIN_CLASS}};
                $value_str = join '; ', @value_array;
            
                $main_class = $value_str;
            }
            if (defined($table_hash{$entry_row}{CATEGORY}))
            {
                @value_array = sort keys %{$table_hash{$entry_row}{CATEGORY}};
                $value_str = join '; ', @value_array;
            
                $category = $value_str;
            }
            
            if ($class_level4 ne '')
            {
                $chosen_name = $class_level4;
            }
            elsif ($sub_class ne '')
            {
                $chosen_name = $sub_class;
            }
            elsif ($main_class ne '')
            {
                $chosen_name = $main_class;
            }
            elsif ($category ne '')
            {
                $chosen_name = $category;
            }
        }
        

        if ($chosen_name ne '')
        {
            $table_hash{$entry_row}{$key}{$chosen_name} = 1;
            $key_count_hash{$key} += 1;
        }
    }
}


@sorted_key_array = sort cmp_keys keys %key_rank_hash;

# print header line
for ($i = 0; $i < $num_keys; $i++)
{
    $key = $sorted_key_array[$i];

    if ($i)
    {
        print "\t";
    }
    
    print "$key";
}
print "\n";


# print table
foreach $entry_row (@entry_row_array)
{
    for ($i = 0; $i < $num_keys; $i++)
    {
        $key       = $sorted_key_array[$i];
        $value_str = '';

        if (defined($table_hash{$entry_row}{$key}))
        {
            @value_array = sort keys %{$table_hash{$entry_row}{$key}};
            
            # warn if we had multiple values stored per entry_row,key pair
            if (@value_array > 1)
            {
                printf STDERR "WARNING -- multiple values for %s:%s\n",
                    $entry_row, $key;
            }
            
            $value_str = join '; ', @value_array;
        }

        if ($i)
        {
            print "\t";
        }

        # replace ; delimiters with |
        $value_str =~ s/\s*;\s+/ \| /g;
        $value_str =~ s/(?<![0-9]);/ \| /g;    # protect lipid abbreviations
    
        print "$value_str";
    }
    print "\n";
}


@array = sort keys %key_with_slash_hash;
foreach $key (@array)
{
    printf STDERR "Substitute / with ; column:\t%s\n", $key;
}
