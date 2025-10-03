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


# sort numeric parts as numbers, not strings
sub cmp_args_alphanumeric
{
    my @array_a = split /([0-9]+)/, $_[0];
    my @array_b = split /([0-9]+)/, $_[1];
    my $count_a = @array_a;
    my $count_b = @array_b;
    my $min_count;
    my $i;
    my $j;
    
    $min_count = $count_a;
    if ($count_b < $min_count)
    {
        $min_count = $count_b;
    }
    
    for ($i = 0; $i < $min_count; $i += 2)
    {
        # even fields sort alphabetically
        if ($array_a[$i] lt $array_b[$i]) { return -1; }
        if ($array_a[$i] gt $array_b[$i]) { return  1; }
        
        # odd fields sort numerically
        $j = $i + 1;
        if ($j < $min_count)
        {
            if ($array_a[$j] < $array_b[$j]) { return -1; }
            if ($array_a[$j] > $array_b[$j]) { return  1; }
        }
    }

    # sort shorter remaining portion first
    if ($count_a < $count_b) { return -1; }
    if ($count_a > $count_b) { return  1; }

    # this shouldn't ever trigger
    return $_[0] cmp $_[1];
}


# handle heavy elements as well
sub cmp_elements
{
    my $ele_a;
    my $ele_b;
    my $heavy_a;
    my $heavy_b;

    $heavy_a = '';
    $heavy_b = '';

    if ($a =~ /^\[([0-9]+)\]/)
    {
        $heavy_a = $1;
    }
    if ($b =~ /^\[([0-9]+)\]/)
    {
        $heavy_b = $1;
    }
    
    $a     =~ /([A-Za-z]+)/;
    $ele_a = $1;
    
    $b     =~ /([A-Za-z]+)/;
    $ele_b = $1;

    # first by element
    if ($ele_a ne $ele_b)
    {
        return $ele_a cmp $ele_b;
    }

    # then put heavy labeled atoms first
    if ($heavy_a ne '' && $heavy_b eq '') { return -1; }
    if ($heavy_b ne '' && $heavy_a eq '') { return  1; }

    # then sort by number of heavy
    if ($heavy_a != $heavy_b)
    {
        return $heavy_a <=> $heavy_b;
    }
   
    return $a cmp $b;
}


# The Hill system specifies C#H#D#, not C#D#H#
#   list all elements in alphabetical order,
#   unless it contains a carbon, then list carbon then hydrogen first
#
# I check for 3-letter elements (all the Uuu's have real symbols by now),
# and elements listed multiple times, and exit early with the original
# formula if such errors are detected.
#
# I don't currently check to see if the given 1- or 2- letter elements are
# valid known elements or not.  I could, but that would require a good bit
# more work than I have time for at the moment.  I'm not *quite* that paranoid
# about the formulas just yet, although part of me still worries about it...
#
sub conform_formula
{
    my $formula_orig = $_[0];
    my $adduct       = $_[1];
    my $formula;
    my $formula_new = '';
    my @match_array;
    my @element_array;
    my @ion_array;    # each +/- group of elements in the adduct
    my $ion;
    my $match;
    my $heavy;
    my $element;
    my $heavy_plus_element;
    my $count;
    my %count_hash = ();
    my %count_adduct_hash = ();
    my $has_carbon_flag = 0;
    my $i;
    
    $formula         = $formula_orig;
    $has_carbon_flag = 0;
    
    # element with number
    @match_array = $formula =~ m/(?:\[[0-9]+\])*[A-Z][a-z]*(?:[0-9]+)*/g;
    foreach $match (@match_array)
    {
        $match =~ /((?:\[[0-9]+\])*)([A-Za-z]+)([0-9]+)*/;

        $heavy   = $1;
        $element = $2;
        $count   = $3;

        if (!defined($heavy))
        {
            $heavy = '';
        }
        if (!defined($count))
        {
            $count = 1;
        }

        if ($element eq 'C')
        {
            $has_carbon_flag = 1;
        }

        $heavy_plus_element = $heavy . $element;

        if (length $element > 2 || defined($count_hash{$heavy_plus_element}))
        {
            printf STDERR "WARNING -- error in formula %s\n", $formula_orig;
            return $formula_orig;
        }
        
        $count_hash{$heavy_plus_element} = $count;
    }
    
    
    # deal with adduct
    if (defined($adduct))
    {
        @ion_array = split /([+-])/, $adduct;
        
        for ($i = 1; $i < @ion_array - 1; $i += 2)
        {
            $add_sub = $ion_array[$i];
            $ion     = $ion_array[$i+1];

            @match_array = $ion =~ m/(?:\[[0-9]+\])*[A-Z][a-z]*(?:[0-9]+)*/g;
            foreach $match (@match_array)
            {
                $match =~ /((?:\[[0-9]+\])*)([A-Za-z]+)([0-9]+)*/;

                $heavy   = $1;
                $element = $2;
                $count   = $3;

                if (!defined($heavy))
                {
                    $heavy = '';
                }
                if (!defined($count))
                {
                    $count = 1;
                }

                # (-) adduct loses a C, so add back into to parent
                if ($add_sub eq '-' && $element eq 'C')
                {
                    $has_carbon_flag = 1;
                }

                $heavy_plus_element = $heavy . $element;

                if (length $element > 2)
                {
                    printf STDERR "WARNING -- error in adduct %s %s\n",
                        $adduct, $ion;

                    return $formula_orig;
                }
                
                if ($add_sub eq '-')
                {
                    $count_hash{$heavy_plus_element} += $count;
                }
                elsif ($add_sub eq '+')
                {
                    $count_hash{$heavy_plus_element} -= $count;
                }
            }
        }
    }
    
    
    @element_array = sort cmp_elements keys %count_hash;
    
    # order all carbons first, followed by hydrogens
    if ($has_carbon_flag)
    {
        # print all carbons first
        foreach $heavy_plus_element (@element_array)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'C')
            {
                $count = $count_hash{$heavy_plus_element};
                
                # skip elements we've subtracted away with the adduct
                if ($count <= 0)
                {
                    next;
                }

                # replace 1 count with blank
                if ($count == 1)
                {
                    $count = '';
                }

                $formula_new .= $heavy_plus_element . $count;
            }
        }
        
        # then all hydrogens that aren't D's
        foreach $heavy_plus_element (@element_array)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'H')
            {
                $count = $count_hash{$heavy_plus_element};

                # skip elements we've subtracted away with the adduct
                if ($count <= 0)
                {
                    next;
                }

                # replace 1 count with blank
                if ($count == 1)
                {
                    $count = '';
                }

                $formula_new .= $heavy_plus_element . $count;
            }
        }
    }

    # order the remaining elements alphabetically, including deuterium
    foreach $heavy_plus_element (@element_array)
    {
        # skip C's and H's we've already placed first
        if ($has_carbon_flag)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'C' || $element eq 'H')
            {
                next;
            }
        }
    
        $count = $count_hash{$heavy_plus_element};

        # skip elements we've subtracted away with the adduct
        if ($count <= 0)
        {
            next;
        }
        
        # replace 1 count with blank
        if ($count == 1)
        {
            $count = '';
        }
        
        $formula_new .= $heavy_plus_element . $count;
    }
    
    return $formula_new;
}


sub check_disallowed_elements
{
    my $formula_orig = $_[0];
    my $formula;
    my $formula_new = '';
    my @match_array;
    my @element_array;
    my @ion_array;    # each +/- group of elements in the adduct
    my $ion;
    my $match;
    my $heavy;
    my $element;
    my $heavy_plus_element;
    my $count;
    my $disallowed_flag;
    my %count_hash = ();
    my $i;
    
    $formula         = $formula_orig;
    $has_carbon_flag = 0;
    
    # element with number
    @match_array = $formula =~ m/(?:\[[0-9]+\])*[A-Z][a-z]*(?:[0-9]+)*/g;
    foreach $match (@match_array)
    {
        $match =~ /((?:\[[0-9]+\])*)([A-Za-z]+)([0-9]+)*/;

        $heavy   = $1;
        $element = $2;
        $count   = $3;

        if (!defined($heavy))
        {
            $heavy = '';
        }
        if (!defined($count))
        {
            $count = 1;
        }

        if ($element eq 'C')
        {
            $has_carbon_flag = 1;
        }

        $heavy_plus_element = $heavy . $element;

        # unsupported formula, skip disfavored check
        if (length $element > 2 || defined($count_hash{$heavy_plus_element}))
        {
            return 0;
        }
        
        $count_hash{$heavy_plus_element} = $count;

        if (!defined($global_count_hash{$heavy_plus_element}))
        {
            $global_count_hash{$heavy_plus_element} = 0;
        }
        $global_count_hash{$heavy_plus_element} += $count;
    }

    @element_array = sort cmp_elements keys %count_hash;

    $disallowed_flag = 0;
    foreach $element (@element_array)
    {
        if ($element =~ /[^SPONCH]/)
        {
            $disallowed_flag = 1;
            last;
        }
    }
    
    return $disallowed_flag;
}


sub cmp_formula
{
    my $formula_1 = $a;
    my $formula_2 = $b;
    my $len1;
    my $len2;
    my $compare;
    
    
    # sort non-polymer formulas first
    if (!($formula_1 =~ /\)[a-z]/) && $formula_2 =~ /\)[a-z]/) { return -1; }
    if (!($formula_2 =~ /\)[a-z]/) && $formula_1 =~ /\)[a-z]/) { return  1; }
    

    # sort longer formulas first (probably have more H's)
    $len1 = length $formula_1;
    $len2 = length $formula_2;
    if ($len1 > $len2) { return -1; }
    if ($len1 < $len2) { return  1; }

    
    # sort larger numbers first
    # InChI can normalize away hydrogens that we want to count
    $compare = cmp_args_alphanumeric($formula_1, $formula_2);
    if ($compare) { return -$compare; }
    
    return 0;
}


sub read_in_darkchem_output
{
    my $filename = $_[0];
    my $adduct   = $_[1];
    my $line;
    my @header_array;
    my @array;
    my $header;
    my $pred_mz_col  = '';
    my $pred_ccs_col = '';
    my $row_col      = '';
    my $pred_mz;
    my $pred_ccs;
    my $row;
    my $i;
    
    open DARKCHEM_OUTPUT, "$filename" or die "ABORT -- cannot open DarkChem output file $filename\n";

    # header line
    $line = <DARKCHEM_OUTPUT>;
    $line =~ s/[\r\n]+$//g;
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        # strip Excel-escaped stuff
        if ($array[$i] =~ /^\"(.*)\"$/)
        {
            $array[$i] = $1;
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\"\"/\"/g;
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
        }
    }
    @header_array = @array;
    
    # search for cols
    for ($i = 0; $i < @array; $i++)
    {
        if ($row_col eq '' && $array[$i] eq 'Row')
        {
            $row_col      = $i;
        }
        elsif ($pred_mz_col eq '' && $array[$i] eq 'prop_000')
        {
            $pred_mz_col  = $i;
        }
        elsif ($pred_ccs_col eq '' && $array[$i] eq 'prop_001')
        {
            $pred_ccs_col = $i;
        }
    }
    
    if ($row_col eq '' || $pred_mz_col eq '' || $pred_ccs_col eq '')
    {
        print "ABORT -- DarkChem output missing Row, prop_000, or prop_001 column(s)\n";
    
        exit(1);
    }
    
    while(defined($line=<DARKCHEM_OUTPUT>))
    {
        $line =~ s/[\r\n]+$//g;
        @array = split /\t/, $line, -1;
        for ($i = 0; $i < @array; $i++)
        {
            # strip Excel-escaped stuff
            if ($array[$i] =~ /^\"(.*)\"$/)
            {
                $array[$i] = $1;
                $array[$i] =~ s/^\s+//;
                $array[$i] =~ s/\s+$//;
                $array[$i] =~ s/\"\"/\"/g;
                $array[$i] =~ s/^\s+//;
                $array[$i] =~ s/\s+$//;
            }
        }
        # pad out any lines that are too short
        for ($i = @array; $i < @header_array; $i++)
        {
            $array[$i] = '';
        }
        
        $row      = $array[$row_col];
        $pred_mz  = $array[$pred_mz_col];
        $pred_ccs = $array[$pred_ccs_col];
        
        # skip rows with missing data
        if ($row eq '' || $pred_mz eq '' || $pred_ccs eq '')
        {
            next;
        }
        
        $darkchem_pred_row_hash{$row}{$adduct}{mz}  = $pred_mz;
        $darkchem_pred_row_hash{$row}{$adduct}{ccs} = $pred_ccs;
    }
    close DARKCHEM_OUTPUT;
}


# begin main()

$darkchem_weights_path =
    '/share/data2/welshea/metabolomics/ccs/darkchem/darkchem-weights';

#$electron_mass      = 0.0005486;
#$h_charged_mass     = 1.0078246 - $electron_mass;   # 1.007276

# all adducts listed here are +1 or -1 charge
$adduct_hash{'[M+H]+'}{offset}       =   1.007276;
$adduct_hash{'[M+Na]+'}{offset}      =  22.989218;
$adduct_hash{'[M-H]-'}{offset}       =  -1.007276;
#$adduct_hash{'[M-H]-'}{requires}{H} =          1;


$temp_darkchem_base   = 'temp_ccs_input_' . $$;
$temp_darkchem_input  = $temp_darkchem_base . '.txt';
$temp_darkchem_output = $temp_darkchem_base . '_darkchem.tsv';


$filename = shift;

open INFILE, $filename or die "can't open input file $filename\n";


# header line
$line = <INFILE>;
$line =~ s/[\r\n]+$//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    # strip Excel-escaped stuff
    if ($array[$i] =~ /^\"(.*)\"$/)
    {
        $array[$i] = $1;
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\"\"/\"/g;
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
}
@header_array = @array;


# take last column encountered instead of first
# since I usually have less-trustworthy --> more-trustworthy left->right
for ($i = 0; $i < @array; $i++)
{
    $field = $array[$i];
    
    # anything with "smiles" in it
    if ($field =~ /SMILES/i)
    {
        $smiles_col = $i;
    }

    # anything with "inchi" in it, but not "key" (skip inchikey, etc.)
    if ($field =~ /InChI/i && !($field =~ /Key/i))
    {
        $inchi_col = $i;
    }
    
    if ($field =~ /formula/i)
    {
        $formula_col = $i;
    }
    
    if ($field =~ /mass/i && $field =~ /mono/i)
    {
        $mass_col = 1;
    }
}


# Normally, we'd prefer InChI over SMILES, because it is more likely to be
# more stereospecific, and is supposed to be unique.  However, several
# HMDB InChI throw OpenBabel errors, so SMILES is cleaner.
#
# Unfortunately, DarkChem trained their models on SMILES strings <= 100,
# so if the input SMILES string is too long, it won't make a prediction.
# Canonical SMILES are often shorter than the HMDB-provided SMILES, so
# we'll let DarkChem do the InChI -> Canonical SMILES conversion for us.
# If you give DarkChem a SMILES string instead of InChI, it does the
# Canonical SMILES conversion and includes it in the output, but we don't get
# the same results as if we had used Open Babel externally to canonicalize
# them.  For whatever reason, we pick up a few extra predictions if we feed
# InChI to DarkChem instead of non-canonical SMILES.  Feeding it Canonical
# SMILES we have converted from the original SMILES strings ourselves with
# Open Babel externally also works, but is slower and more tedious to program
# around, so we'll just let Dark Chem convert from InChI.
#
# Note to future self:  Inchified SMILES are better than Universal SMILES.
# Open Babel Canonical SMILES is *NOT* the same as Inchified SMILES, and
# does not work quite as well (but is still better than Universal SMILES).
# Open Babel Canonical SMILES applied to InChI *might* be identical to
# Inchified Smiles ??
#
# See O'Boyle 2012:
#   https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-22

$smiles_or_inchi = '';
if (defined($inchi_col))
{
    $smiles_or_inchi = 'inchi';
    $struct_col      = $inchi_col;
}
elsif (defined($smiles_col))
{
    $smiles_or_inchi = 'smiles';
    $struct_col      = $smiles_col;
}


# read in data file
$row = 0;
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+$//g;
    @array = split /\t/, $line, -1;
    for ($i = 0; $i < @array; $i++)
    {
        # strip Excel-escaped stuff
        if ($array[$i] =~ /^\"(.*)\"$/)
        {
            $array[$i] = $1;
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\"\"/\"/g;
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
        }
    }
    # pad out any lines that are too short
    for ($i = @array; $i < @header_array; $i++)
    {
        $array[$i] = '';
    }

    # replace line with cleaned line
    $line = join "\t", @array;
    
    $smiles        = '';
    $inchi         = '';
    $formula_orig  = '';
    $formula       = '';
    $formula_inchi = '';
    $mass          = '';
    
    if (defined($smiles_col) && $array[$smiles_col] =~ /[A-Za-z]/)
    {
        $smiles = $array[$smiles_col];
    }
    if (defined($inchi_col) && $array[$inchi_col] =~ /[A-Za-z]/)
    {
        $inchi = $array[$inchi_col];
    }
    if (defined($formula_col) && $array[$formula_col] =~ /[A-Za-z]/)
    {
        $formula_orig = $array[$formula_col];
    }
    if (defined($mass_col) && is_number($array[$mass_col]))
    {
        $mass = $array[$mass_col];
    }


    $formula_inchi = '';
    if ($inchi ne '' && $inchi =~ /\/([^\/]+)\//)
    {
        $formula_inchi = $1;
    }


    # skip anything with ionic bonds
    if ($smiles ne '' && $smiles =~ /\./)
    {
        next;
    }
    if ($formula_orig ne '' && $formula_orig =~ /\./)
    {
        next;
    }
    if ($formula_inchi ne '' && $formula_inchi =~ /\./)
    {
        next;
    }
    
    
    # skip anything with missing structure in chosen column
    if (!($array[$struct_col] =~ /[A-Za-z]/))
    {
        next;
    }
    
    
    # conform formula, don't conform polymers (not supported by my code)
    $formula = $formula_orig;
    if ($formula_orig ne '' && !($formula_orig =~ /\)[a-z]/))
    {
        $formula = conform_formula($formula_orig);
    }
    if ($formula_inchi ne '' && !($formula_inchi =~ /\)[a-z]/))
    {
        $formula_inchi = conform_formula($formula_inchi);
    }
    
    
    @formula_array = ();
    $num_formulas = 0;
    if ($formula ne '')
    {
        $formula_array[$num_formulas++] = $formula;
    }
    if ($formula_inchi ne '')
    {
        $formula_array[$num_formulas++] = $formula_inchi;
    }

    $formula_chosen = '';
    if ($num_formulas)
    {
        @formula_array = sort cmp_formula @formula_array;
        $formula_chosen = $formula_array[0];
    }


    # skip formulas with disallowed elements (only allow SPONCH)
    # also skip unsupported formulas syntax, such as polymers
    # NOTE -- we've already chosen between polymer and non-polymer formulas
    if ($formula_chosen ne '' && check_disallowed_elements($formula_chosen))
    {
        next;
    }
    
    
    # skip anything without a hydrogen
    # these are usually industrial chemicals,
    # and unlikely to ever be relevant to our studies
    #
    # Also, most CCS predictors fail to recognize that molecules without H
    # cannot be deprotonated, and thus happily predict CCS values for
    # impossible [M-H]- adducts...
    #
    if (1 && $formula_chosen ne '' && !($formula_chosen =~ /H/))
    {
        next;
    }
    

    # sanity check formula and inchi
    # skip polymers
    if ($formula ne '' && $formula_inchi && !($formula =~ /\)[a-z]/))
    {
        # strip hydrogens, since HMDB is often off by some hydrogens
        $formula_no_h       =  $formula;
        $formula_inchi_no_h =  $formula_inchi;
        $formula_no_h       =~ s/[H][0-9]*//g;
        $formula_inchi_no_h =~ s/[H][0-9]*//g;
        
        if ($formula_no_h ne $formula_inchi_no_h)
        {
            printf STDERR "WARNING -- %s\t%s\t%s\n",
                $array[0], $formula, $formula_inchi;
        }
    }


    $data_array[$row]{id}      = $array[0];
    $data_array[$row]{smiles}  = $smiles;
    $data_array[$row]{inchi}   = $inchi;
    $data_array[$row]{formula} = $formula_chosen;
    $data_array[$row]{mass}    = $mass;
    $data_array[$row]{line}    = $line;


    if (0)
    {
        if ($num_formulas == 2 &&
            $formula_chosen ne '' &&
            $formula_array[0] ne '' &&
            $formula_array[1] ne '' &&
            $formula_array[0] ne $formula_array[1])
        {
            printf STDERR "%s\t%s\t%s\t!=\t%s\n",
                $array[0], $row,
                $formula_array[0], $formula_array[1];
        }
    }

    
    $row++;
}

$num_rows = $row;



# generate input file for DarkChem
open DARKCHEM_INPUT, ">$temp_darkchem_input" or die "ABORT -- cannot open file $temp_darkchem_input\n";
printf DARKCHEM_INPUT "%s\t%s\t%s\t%s",
    'Identifier', 'Row', 'Mass', 'Formula';

if ($smiles_or_inchi eq 'inchi')
{
    print DARKCHEM_INPUT "\tInChI";
}
else
{
    print DARKCHEM_INPUT "\tSMILES";
}
print DARKCHEM_INPUT "\n";

for ($row = 0; $row < $num_rows; $row++)
{
    # we've already skipped any rows with missing structures earlier

    $id      = $data_array[$row]{id};
    $mass    = $data_array[$row]{mass};
    $formula = $data_array[$row]{formula};

    print DARKCHEM_INPUT "$id\t$row\t$mass\t$formula";

    $struct = '';
    if ($smiles_or_inchi eq 'inchi')
    {
        $struct = $data_array[$row]{inchi};
    }
    else
    {
        $struct = $data_array[$row]{smiles};
    }
    
    print DARKCHEM_INPUT "\t$struct";
    print DARKCHEM_INPUT "\n";
}
close DARKCHEM_INPUT;



# run DarkChem on M+H
$cmd_str = sprintf "darkchem predict prop \"%s\" \"%s/%s\" 2>&1 > /dev/null",
           $temp_darkchem_input,
           $darkchem_weights_path,
           'protonated';
`$cmd_str`;
read_in_darkchem_output($temp_darkchem_output, '[M+H]+');


# run DarkChem on M+Na
$cmd_str = sprintf "darkchem predict prop \"%s\" \"%s/%s\" 2>&1 > /dev/null",
           $temp_darkchem_input,
           $darkchem_weights_path,
           'sodiated';
`$cmd_str`;
read_in_darkchem_output($temp_darkchem_output, '[M+Na]+');


# run DarkChem on M-H
$cmd_str = sprintf "darkchem predict prop \"%s\" \"%s/%s\" 2>&1 > /dev/null",
           $temp_darkchem_input,
           $darkchem_weights_path,
           'deprotonated/';
`$cmd_str`;
read_in_darkchem_output($temp_darkchem_output, '[M-H]-');



$adduct_hash{'[M+H]+'} = 1;
$adduct_hash{'[M+Na]+'} = 1;
$adduct_hash{'[M-H]-'} = 1;

@adduct_array = sort keys %adduct_hash;



# print header
$header_line  = join "\t", @header_array;
$header_line .= "\tformula_chosen";
foreach $adduct (@adduct_array)
{
    $header_line .= sprintf "\tCCS: %s", $adduct;
}
if (0)
{
    foreach $adduct (@adduct_array)
    {
        $header_line .= sprintf "\tDarkChem m/z: %s", $adduct;
    }
}
print "$header_line\n";



# output table of CCS predictions
for ($row = 0; $row < $num_rows; $row++)
{
    $formula_chosen = $data_array[$row]{formula};

    $line = $data_array[$row]{line};
    $line .= "\t$formula_chosen";

    foreach $adduct (@adduct_array)
    {
        $ccs = $darkchem_pred_row_hash{$row}{$adduct}{ccs};
        if (!defined($ccs))
        {
            $ccs = '';
        }
        
        # blank out impossible [M-H]- adducts
        # DarkChem, DeepCCS, AllCCS all fail on this
        $formula = $data_array[$row]{formula};
        if ($adduct eq '[M-H]-' && $formula ne '' && !($formula =~ /H/))
        {
            $ccs = '';
        }
    
        $line .= sprintf "\t%s", $ccs;
    }
    if (0)
    {
      foreach $adduct (@adduct_array)
      {
        $mz = $darkchem_pred_row_hash{$row}{$adduct}{mz};
        if (!defined($mz))
        {
            $mz = '';
        }

        # blank out impossible [M-H]- adducts
        # DarkChem, DeepCCS, AllCCS all fail on this
        $formula = $data_array[$row]{formula};
        if ($adduct eq '[M-H]-' && $formula ne '' && !($formula =~ /H/))
        {
            $mz = '';
        }
    
        $line .= sprintf "\t%s", $mz;
      }
    }

    print "$line\n";
}


# remove temp files
`rm $temp_darkchem_input`;
`rm $temp_darkchem_output`;
