#!/usr/bin/perl -w

use Scalar::Util qw(looks_like_number);
use POSIX;
use File::Basename;


sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not.
        #
        # Perl treats NaN or NaNs, and various mixed caps, as numbers.
        # Weird that not-a-number is a number... but it is so that
        # it can do things like nan + 1 = nan, so I guess it makes sense
        #
        if ($_[0] =~ /^[-+]*(Inf|NaN)/i)
        {
            return 0;
        }
        
        return 1;
    }

    # optional + or - sign at beginning
    # then require either:
    #  a number followed by optional comma stuff, then optional decimal stuff
    #  mandatory decimal, followed by optional digits
    # then optional exponent stuff
    #
    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    if ($_[0] =~ /^([-+]?)([0-9]+(,[0-9]{3,})*\.?[0-9]*|\.[0-9]*)([Ee]([-+]?[0-9]+))?$/)
    {
        # current REGEX can treat '.' as a number, check for that
        if ($_[0] eq '.')
        {
            return 0;
        }
        
        return 1;
    }
    
    return 0;
}


# requires tab-delimited text input, must convert semi-colon to tab beforehand
#
# SCILS fails to escape ; present in the Name field,
# resulting in corrupt .csv exports
#
# For now, require NeutralMass to immediately follow Name column,
# since we use that to check when things have gone wrong

$infile = shift;	# already converted to dab-delimited

if (!defined($infile))
{
    $infile = '-';
}

open MALDI_CCS, "$infile" or die "ABORT -- cannot open file $infile\n";
$row = -1;


# read in MALDI ccs file
# jump down to header line
while(defined($line=<MALDI_CCS>))
{
    $row++;

    # pass comment line through as-is
    $c = substr $line, 0, 1;
    if ($c eq '#')
    {
        print $line;
        next;
    }

    # found likely header line
    if ($c =~ /\S/)
    {
        print $line;
        last;
    }
}

$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    if ($array[$i] =~ s/^\"(.*)\"$/$1/)
    {
        $array[$i] =~ s/\"\"/\"/g;
    }
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    
    $field = $array[$i];

    if ($field =~ /^Name$/i)
    {
        $name_col = $i;
    }
}
#$maldi_ccs_header_line = join "\t", @array;
#$num_header_cols = @array;

$numeric_col    = $name_col + 1;
$numeric_header = $array[$numeric_col];

$has_required_col_order_flag = 0;
if (defined($numeric_header) &&
    ($numeric_header =~ /Mass/i ||
     $numeric_header =~ /Intensity/i))
{
    $has_required_col_order_flag = 1;
}


if ($has_required_col_order_flag == 0)
{
    printf STDERR "WARNING -- Name followed by %s col, cannot uncorrupt\n";
}


# read in the MALDI ccs data
while(defined($line=<MALDI_CCS>))
{
    $row++;

    # unsupported column order, pass though as-is
    if ($has_required_col_order_flag == 0)
    {
        print $line;
        next;
    }

    # pass comment line through as-is
    $c = substr $line, 0, 1;
    if ($c eq '#')
    {
        print $line;
        next;
    }

    @array = split /\t/, $line, -1;

    #for ($i = 0; $i < @array; $i++)
    #{
    #    $array[$i] =~ s/^\s+//;
    #    $array[$i] =~ s/\s+$//;
    #    if ($array[$i] =~ s/^\"(.*)\"$/$1/)
    #    {
    #        $array[$i] =~ s/\"\"/\"/g;
    #    }
    #    $array[$i] =~ s/^\s+//;
    #    $array[$i] =~ s/\s+$//;
    #}
    
    $numeric = $array[$numeric_col];
    $name    = $array[$name_col];
    
    # uh oh, fields do not exist
    if (!defined($numeric) || !defined($name))
    {
        printf STDERR "WARNING -- row %d too short, passing through unchanged\n";
    }
    
    # NeutralMass should always be present, and always be numeric
    #
    # Corrupted Name fields from unescaped embedded ; never result in
    # purely numeric data (so far).
    #
    # So, we'll use this to detect and correct corrupted fields.
    
    # everything is OK, pass through as-is
    # NeutralMass is a number, or Name is empty or all whitespace (blank)
    if (is_number($numeric) || !($name =~ /\S/))
    {
        print $line;
        next;
    }
    
    # concatenate corrupted Name fields until we reach a numeric field
    $name_str     = $name;
    $col_scan_end = $name_col + 1;
    for ($i = $name_col + 1; $i < @array; $i++)
    {
        $field = $array[$i];
        
        if (!defined($field))
        {
            printf STDERR "WARNING -- row %d too short, row may be truncated\n";

            $col_scan_end = $i - 1;
            last;
        }
    
        # found what we assume is the NeutralMass value
        if (is_number($field))
        {
            $col_scan_end = $i;
            last;
        }
        
        $name_str .= ';' . $field;
    }
    
    if ($col_scan_end == $name_col + 1)
    {
        printf STDERR "WARNING -- row %d reached code it shouldn't have...\n";

        print $line;
        next;
    }
    

    printf STDERR "Uncorrupting row %d:\t%s\n",
        $row, $name_str;


    # store up until Name as-is
    @array_new = ();
    for ($i = 0; $i < $name_col; $i++)
    {
        $array_new[$i] = $array[$i];
    }

    # store the new Name
    $array_new[$i++] = $name_str;

    # shift everything beyond Name
    for ($j = $col_scan_end; $j < @array; $j++)
    {
        $field = $array[$j];
        
        if (!defined($field))
        {
            printf STDERR "WARNING -- row %d too short, row may be truncated\n";
            last;
        }
        
        $array_new[$i++] = $field;
    }

    $line_new = join "\t", @array_new;

    print $line_new;
}
close MALDI_CCS;
