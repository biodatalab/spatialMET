#!/usr/bin/perl -w


# 10-20 PPM for m/z, 3% for CCS

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


# lookup table for converting numbers into ordered greek letters
$number_letter_hash{1} = 'a';
$number_letter_hash{2} = 'b';
$number_letter_hash{3} = 'g';
$number_letter_hash{4} = 'd';
$number_letter_hash{5} = 'e';
$number_letter_hash{6} = 'z';
$number_letter_hash{7} = 'h';

# common Greek and punctuation unicode seen in metabolite names
# convert them to their nearest ASCII equivalent
#
# Greek
$unicode_to_ascii_hash{"\x{0391}"} = 'A';
$unicode_to_ascii_hash{"\x{0392}"} = 'B';
$unicode_to_ascii_hash{"\x{0393}"} = 'G';
$unicode_to_ascii_hash{"\x{0394}"} = 'D';
$unicode_to_ascii_hash{"\x{0395}"} = 'E';
$unicode_to_ascii_hash{"\x{0396}"} = 'Z';
$unicode_to_ascii_hash{"\x{0397}"} = 'H';
$unicode_to_ascii_hash{"\x{039B}"} = 'L';
$unicode_to_ascii_hash{"\x{03A6}"} = 'Phi';
$unicode_to_ascii_hash{"\x{03A8}"} = 'Psi';
$unicode_to_ascii_hash{"\x{03A9}"} = 'O';
$unicode_to_ascii_hash{"\x{03B1}"} = 'a';
$unicode_to_ascii_hash{"\x{03B2}"} = 'b';
$unicode_to_ascii_hash{"\x{03B3}"} = 'g';
$unicode_to_ascii_hash{"\x{03B4}"} = 'd';
$unicode_to_ascii_hash{"\x{03B5}"} = 'e';
$unicode_to_ascii_hash{"\x{03B6}"} = 'z';
$unicode_to_ascii_hash{"\x{03B7}"} = 'h';
$unicode_to_ascii_hash{"\x{03BB}"} = 'l';
$unicode_to_ascii_hash{"\x{03C6}"} = 'phi';
$unicode_to_ascii_hash{"\x{03C8}"} = 'psi';
$unicode_to_ascii_hash{"\x{03C9}"} = 'o';


# Latin-1 supplement letters block
if (0)
{
$unicode_to_ascii_hash{"\x{00C0}"} = 'A';
$unicode_to_ascii_hash{"\x{00C1}"} = 'A';
$unicode_to_ascii_hash{"\x{00C2}"} = 'A';
$unicode_to_ascii_hash{"\x{00C3}"} = 'A';
$unicode_to_ascii_hash{"\x{00C4}"} = 'A';
$unicode_to_ascii_hash{"\x{00C5}"} = 'A';
$unicode_to_ascii_hash{"\x{00C6}"} = 'AE';
$unicode_to_ascii_hash{"\x{00C7}"} = 'C';
$unicode_to_ascii_hash{"\x{00C8}"} = 'E';
$unicode_to_ascii_hash{"\x{00C9}"} = 'E';
$unicode_to_ascii_hash{"\x{00CA}"} = 'E';
$unicode_to_ascii_hash{"\x{00CB}"} = 'E';
$unicode_to_ascii_hash{"\x{00CC}"} = 'I';
$unicode_to_ascii_hash{"\x{00CD}"} = 'I';
$unicode_to_ascii_hash{"\x{00CE}"} = 'I';
$unicode_to_ascii_hash{"\x{00CF}"} = 'I';
$unicode_to_ascii_hash{"\x{00D0}"} = 'D';
$unicode_to_ascii_hash{"\x{00D1}"} = 'N';
$unicode_to_ascii_hash{"\x{00D2}"} = 'O';
$unicode_to_ascii_hash{"\x{00D3}"} = 'O';
$unicode_to_ascii_hash{"\x{00D4}"} = 'O';
$unicode_to_ascii_hash{"\x{00D5}"} = 'O';
$unicode_to_ascii_hash{"\x{00D6}"} = 'O';
$unicode_to_ascii_hash{"\x{00D7}"} = 'x';    # multiplication
$unicode_to_ascii_hash{"\x{00D8}"} = 'O';
$unicode_to_ascii_hash{"\x{00D9}"} = 'U';
$unicode_to_ascii_hash{"\x{00DA}"} = 'U';
$unicode_to_ascii_hash{"\x{00DB}"} = 'U';
$unicode_to_ascii_hash{"\x{00DC}"} = 'U';
$unicode_to_ascii_hash{"\x{00DD}"} = 'Y';
$unicode_to_ascii_hash{"\x{00DE}"} = 'TH';
$unicode_to_ascii_hash{"\x{00DF}"} = 'ss';
$unicode_to_ascii_hash{"\x{00E0}"} = 'a';
$unicode_to_ascii_hash{"\x{00E1}"} = 'a';
$unicode_to_ascii_hash{"\x{00E2}"} = 'a';
$unicode_to_ascii_hash{"\x{00E3}"} = 'a';
$unicode_to_ascii_hash{"\x{00E4}"} = 'a';
$unicode_to_ascii_hash{"\x{00E5}"} = 'a';
$unicode_to_ascii_hash{"\x{00E6}"} = 'ae';
$unicode_to_ascii_hash{"\x{00E7}"} = 'c';
$unicode_to_ascii_hash{"\x{00E8}"} = 'e';
$unicode_to_ascii_hash{"\x{00E9}"} = 'e';
$unicode_to_ascii_hash{"\x{00EA}"} = 'e';
$unicode_to_ascii_hash{"\x{00EB}"} = 'e';
$unicode_to_ascii_hash{"\x{00EC}"} = 'i';
$unicode_to_ascii_hash{"\x{00ED}"} = 'i';
$unicode_to_ascii_hash{"\x{00EE}"} = 'i';
$unicode_to_ascii_hash{"\x{00EF}"} = 'i';
$unicode_to_ascii_hash{"\x{00F0}"} = 'd';
$unicode_to_ascii_hash{"\x{00F1}"} = 'n';
$unicode_to_ascii_hash{"\x{00F2}"} = 'o';
$unicode_to_ascii_hash{"\x{00F3}"} = 'o';
$unicode_to_ascii_hash{"\x{00F4}"} = 'o';
$unicode_to_ascii_hash{"\x{00F5}"} = 'o';
$unicode_to_ascii_hash{"\x{00F6}"} = 'o';
$unicode_to_ascii_hash{"\x{00F7}"} = '/';    # division
$unicode_to_ascii_hash{"\x{00F8}"} = 'o';
$unicode_to_ascii_hash{"\x{00F9}"} = 'u';
$unicode_to_ascii_hash{"\x{00FA}"} = 'u';
$unicode_to_ascii_hash{"\x{00FB}"} = 'u';
$unicode_to_ascii_hash{"\x{00FC}"} = 'u';
$unicode_to_ascii_hash{"\x{00FD}"} = 'y';
$unicode_to_ascii_hash{"\x{00FE}"} = 'th';
$unicode_to_ascii_hash{"\x{00FF}"} = 'y';
# OE from Latin Extended-A
$unicode_to_ascii_hash{"\x{0152}"} = 'OE';
$unicode_to_ascii_hash{"\x{0153}"} = 'oe';
}
# only the most common ones, from synonyms field
# name field can contain junk, which I may want to try to salvage later
else
{
$unicode_to_ascii_hash{"\x{00E4}"} = 'a';
$unicode_to_ascii_hash{"\x{00E8}"} = 'e';
$unicode_to_ascii_hash{"\x{00E9}"} = 'e';
$unicode_to_ascii_hash{"\x{00EF}"} = 'i';
$unicode_to_ascii_hash{"\x{00F4}"} = 'o';
$unicode_to_ascii_hash{"\x{00F6}"} = 'o';
$unicode_to_ascii_hash{"\x{00FC}"} = 'u';
}


# +/-, dashes, and quotes
$unicode_to_ascii_hash{"\x{00B1}"} = '+-';   # synonyms use (+-), not (+/-)
$unicode_to_ascii_hash{"\x{2010}"} = '-';    # name only, real (50 of them)
$unicode_to_ascii_hash{"\x{2011}"} = '-';
$unicode_to_ascii_hash{"\x{2012}"} = '-';
$unicode_to_ascii_hash{"\x{2013}"} = '-';
$unicode_to_ascii_hash{"\x{2014}"} = '-';
$unicode_to_ascii_hash{"\x{2015}"} = '-';
$unicode_to_ascii_hash{"\x{2192}"} = '-';
$unicode_to_ascii_hash{"\x{2212}"} = '-';
#$unicode_to_ascii_hash{"\x{00A8}"} = '"';    # name only, corrupt text
$unicode_to_ascii_hash{"\x{02B9}"} = "'";
$unicode_to_ascii_hash{"\x{02BA}"} = '"';
$unicode_to_ascii_hash{"\x{2018}"} = "'";
$unicode_to_ascii_hash{"\x{2019}"} = "'";
#$unicode_to_ascii_hash{"\x{201A}"} = "'";    # name only, corrupt text
$unicode_to_ascii_hash{"\x{201B}"} = "'";
$unicode_to_ascii_hash{"\x{201C}"} = '"';
$unicode_to_ascii_hash{"\x{201D}"} = '"';
$unicode_to_ascii_hash{"\x{201E}"} = '"';
$unicode_to_ascii_hash{"\x{201F}"} = '"';
$unicode_to_ascii_hash{"\x{2032}"} = "'";
$unicode_to_ascii_hash{"\x{2033}"} = '"';
$unicode_to_ascii_hash{"\x{2035}"} = "'";
$unicode_to_ascii_hash{"\x{2034}"} = "'''";
$unicode_to_ascii_hash{"\x{2036}"} = '"';
$unicode_to_ascii_hash{"\x{2037}"} = "'''";   # triple prime
$unicode_to_ascii_hash{"\x{2057}"} = '""';    # quadruple prime
$unicode_to_ascii_hash{"\x{301D}"} = '"';
$unicode_to_ascii_hash{"\x{301E}"} = '"';
$unicode_to_ascii_hash{"\x{301F}"} = '"';

# other punctuation
# should be removed, shouldn't be there in the first place
$unicode_to_ascii_hash{"\x{00AB}"} = '';   # <<;  alpha,<<gamma>>-butadiene
$unicode_to_ascii_hash{"\x{00BB}"} = '';   # >>;  alpha,<<gamma>>-butadiene
$unicode_to_ascii_hash{"\x{2020}"} = '';   # dagger;        end of HMDB0240697
$unicode_to_ascii_hash{"\x{2021}"} = '';   # double dagger; end of HMDB0240697


# superscript numbers
$unicode_to_ascii_hash{"\x{00B2}"} = '2';
$unicode_to_ascii_hash{"\x{00B3}"} = '3';
$unicode_to_ascii_hash{"\x{00B9}"} = '1';
$unicode_to_ascii_hash{"\x{2070}"} = '0';
$unicode_to_ascii_hash{"\x{2074}"} = '4';
$unicode_to_ascii_hash{"\x{2075}"} = '5';
$unicode_to_ascii_hash{"\x{2076}"} = '6';
$unicode_to_ascii_hash{"\x{2077}"} = '7';
$unicode_to_ascii_hash{"\x{2078}"} = '8';
$unicode_to_ascii_hash{"\x{2079}"} = '9';
# subscript numbers
$unicode_to_ascii_hash{"\x{2080}"} = '0';
$unicode_to_ascii_hash{"\x{2081}"} = '1';
$unicode_to_ascii_hash{"\x{2082}"} = '2';
$unicode_to_ascii_hash{"\x{2083}"} = '3';
$unicode_to_ascii_hash{"\x{2084}"} = '4';
$unicode_to_ascii_hash{"\x{2085}"} = '5';
$unicode_to_ascii_hash{"\x{2086}"} = '6';
$unicode_to_ascii_hash{"\x{2087}"} = '7';
$unicode_to_ascii_hash{"\x{2088}"} = '8';
$unicode_to_ascii_hash{"\x{2089}"} = '9';
# zero-width spaces, should be removed
$unicode_to_ascii_hash{"\x{200B}"} = '';
$unicode_to_ascii_hash{"\x{FEFF}"} = '';    # name only, Bosutinib

# U+00A0 (non-breaking space)
#
# HMDB0062476
#    GalNAc(3S)-GlcA-Gal-Gal-Xyl??
# Non-breaking space at the end
# This occurs several times, likely a copy/paste error
$unicode_to_ascii_hash{"\x{00A0}"} = '';

# U+00AC (NOT symbol)
#
# appears to be inserted junk in front of +/-
# example: HMDB0303381 (+/-)-Isobornyl acetate
#          https://foodb.ca/compounds/FDB012445
# Or part of corrupted multibyte unicode:
#     HMDB0251069
#     HMDB0250632

# bogus character, probably supposed to be alpha?
# it only occurs once in all of the HMDB compound name fields
# HMDB0242122 ?-D-galactopyranoside, ethyl
$unicode_to_ascii_hash{"\x{FFFD}"} = 'a';


sub unicode_to_ascii
{
    my $value = $_[0];
    my $len;
    my $string_new;
    
    if ($value =~ /[\x80-\xFF]/)
    {
        # first, decode the unicode string
        # into single characters, so substr works correctly
        utf8::decode($value);

        #if ($value =~ /[\x{0370}-\x{03ff}]/)
        #{
        #    $temp = $value;
        #    utf8::encode($temp);
        #    printf STDERR "$accession\t$temp\n";
        #}

        ## HACK -- {NOT}+/-
        $value =~ s/\(\x{00AC}\x{00B1}\)/\(\x{00B1}\)/g;
        
        ## HACK -- HMDB0304547
        ## corrupted omega
        ## thank you python ftfy package for confirming!
        $value =~ s/\x{0153}\x{00E2}/\x{03C9}/g;
        
        ## HACK -- HMDB0304570 HMDB0304569
        ## appears to be corrupted '
        ## thank you python ftfy package!
        $value =~ s/\x{201A}\x{00C4}\x{2264}/\'/g;
        
        ## HMDB0300900    (2E_4Z)???\decadienoyl-CoA
        ##                http://qpmf.rx.umaryland.edu/PAMDB?MetID=PAMDB001410
        ##                should probably be (2E_4Z)-decadienoyl-CoA
        ## remove offending unicode entirely, plus the following backslash
        $value =~ s/\x{00D4}\x{00F8}\x{03A9}\\//g;

        ## HMDB0251069    2,2???-(Hydroxynitrosohydrazino)bis-ethanamine
        ##                bloodexposome.org: 2,2'-(Hydroxynitrosohydrazino)...
        ## remove offending unicode entirely
        $value =~ s/\x{201A}\x{00C4}\x{00F6}\x{221A}\x{00D1}\x{221A}\x{2202}\?//g;
        $value =~ s/\x{201A}\x{00C4}\x{00F6}\x{221A}\x{00A2}\x{00AC}\x{00DF}//g;
        
        ## HMDB0250632    ... cyclic (3?5)-disulfide
        ##                bloodexposome.org: ... cyclic (35)-disulfide
        ## remove offending unicode entirely
        $value =~ s/\x{00AC}\x{00A8}\x{00AC}\x{00AE}//g;
        $value =~ s/\x{201A}\x{00E0}\x{00F6}\x{221A}\x{00FA}//g;


        $string_new = '';
        $len        = length $value;
        for ($j = 0; $j < $len; $j++)
        {
            $c_single = substr $value, $j, 1;
            $c_new = $unicode_to_ascii_hash{$c_single};

            if (defined($c_new))
            {
                $string_new .= $c_new;
            }
            else
            {
                $string_new .= $c_single;
            }
        }
        $value = $string_new;
        
        # remove leading/trailing whitespace
        $value =~ s/^\s+//;
        $value =~ s/\s+$//;

        # encode it back again, since input is multi-byte chars
        utf8::encode($value);
    }
    
    return $value;
}

# generate conformed name (strip punctuation, capitalization, etc.)
sub conform_name
{
    my $name       = $_[0];
    my $name_len;
    my $half_name_len;
    my $name_half1 = '';
    my $name_half2 = '';

    my @temp_array;
    my $i;
    
    # lowercase everything
    $name = lc $name;
    
    # convert unicode to ASCII
    $name = unicode_to_ascii($name);
    
    # convert _ to space, so \w and \b won't match on them
    $name =~ s/_/ /g;
    
    # remove likely abbreviations or synonyms from end
    $name =~ s/ +\([^0-9()]+\)$//;

    # replace Greek letters at non-letter boundaries
    # example: HMDB0000708 Glycoursodeoxycholic acid
    $name =~ s/(?<![A-Za-z])alpha(?![A-Za-z])/a/g;
    $name =~ s/(?<![A-Za-z])beta(?![A-Za-z])/b/g;
    $name =~ s/(?<![A-Za-z])gamma(?![A-Za-z])/g/g;
    $name =~ s/(?<![A-Za-z])delta(?![A-Za-z])/d/g;
    $name =~ s/(?<![A-Za-z])epsilon(?![A-Za-z])/e/g;
    $name =~ s/(?<![A-Za-z])zeta(?![A-Za-z])/z/g;
    $name =~ s/(?<![A-Za-z])eta(?![A-Za-z])/h/g;
    
    # replace single numbers with romanized greek letters
    #    2-aminoethylphosphonate --> b-aminoethylphosphonate
    @temp_array = split /\b([1-7])\b/, $name;
    for ($i = 1; $i < @temp_array; $i += 2)
    {
        $temp_array[$i] = $number_letter_hash{$temp_array[$i]};
    }
    $name = join '', @temp_array;

    # conform acids
    # example of "acids": HMDB0001202
    #
    $name =~ s/\bacids\b/acid/g;
    $name =~ s/acids*,(.*)ic$/$1ic acid/; # reorder weird MeSH, HMDB notation
    $name =~ s/anoic/yric/g;              # Butanoic acid --> Butyric acid
    $name =~ s/anoate\b/yrate/g;          # Butanoate     --> Butyrate
    $name =~ s/ic acids*\b/ate/g;         # Glutamic acid --> Glutamate
    $name =~ s/ates\b/ate/g;              # Benzoates     --> Benzoate
    
    # strip the L- from L-aminoacids
    # only on word boundary, so we don't strip DL-aminoacid
    $name =~ s/\bl-(alanine)/$1/g;
    $name =~ s/\bl-(arginine)/$1/g;
    $name =~ s/\bl-(asparagine)/$1/g;
    $name =~ s/\bl-(aspartic acid)/$1/g;
    $name =~ s/\bl-(cysteine)/$1/g;
    $name =~ s/\bl-(glutamine)/$1/g;
    $name =~ s/\bl-(glutamic acid)/$1/g;
    $name =~ s/\bl-(glycine)/$1/g;    # L-Glycine is used, but not by itself
    $name =~ s/\bl-(histidine)/$1/g;
    $name =~ s/\bl-(isoleucine)/$1/g;
    $name =~ s/\bl-(leucine)/$1/g;
    $name =~ s/\bl-(lysine)/$1/g;
    $name =~ s/\bl-(methionine)/$1/g;
    $name =~ s/\bl-(phenylalanine)/$1/g;
    $name =~ s/\bl-(proline)/$1/g;
    $name =~ s/\bl-(serine)/$1/g;
    $name =~ s/\bl-(threonine)/$1/g;
    $name =~ s/\bl-(tryptophan)/$1/g;          # also Tryptophanamide
    $name =~ s/\bl-(tyrosine)/$1/g;
    $name =~ s/\bl-(valine)/$1/g;

    # sulfid/sulfide/sulphid/sulphide
    # HMDB0042033 Thiodiglycol is the only entry with sulfid/sulphid
    # so, replace sulfid/sulphid with sulfide, since sulfid/sulphid is odd
    # sulfide is kept over sulphide due to fewer characters
    #
    $name =~ s/sulph/sulf/g;
    $name =~ s/sulfid\b/sulfide/g;
    #$name =~ s/sulphid\b/sulfide/g;
    #$name =~ s/sulphide\b/sulfide/g;

    # "ic acid" / "ate"
    $name =~ s/\bl-(aspartate)\b/$1/g;
    $name =~ s/\bl-(glutamate)\b/$1/g;

    # artificial, non-Human amino acids or dipeptides
    $name =~ s/\bl-(cysteinylglycine)/$1/g;    # Cysteinylglycine
    $name =~ s/\bl-(homocysteine)/$1/g;        # Homocysteine
    $name =~ s/\bl-(norleucine)/$1/g;          # Norleucine
    $name =~ s/\bl-(selenomethionine)/$1/g;    # Selenomethionine
    $name =~ s/\bl-(anserine)/$1/g;            # Anserine
    $name =~ s/\bl-(homoserine)/$1/g;          # Homoserine
    $name =~ s/\bl-(allothreonine)/$1/g;       # Allothreonine
    $name =~ s/\bl-(norvaline)/$1/g;           # Norvaline
    
    # ethyl, methyl, etc.
    $name =~ s/thane/thyl/g;          # (2-Aminoethane)phosphonic acid
                                      # 2-Aminoethylphosphonate

    # mono is redundant and often left out in synonyms
    #
    # WARNING:  conform to the same string, but methyl in different places
    #   HMDB0000752  Methylglutaric acid
    #   HMDB0000858  Monomethyl glutaric acid
    #
    $name =~ s/mono(\S)/$1/g;
    
    # replace (+/-) with (+-)
    $name =~ s/\(\+\/\-\)/\(\+\-\)/g;

    # keep only D,L,DL when together with (+),(-),(+-)
    #
    # examples: D-(-)-Arabinose
    #           D-(+)-Galactosamine
    #           D-(+)-Galacturonic acid
    #           D-(+)-Glucosamine
    #
    #swap them around if +/- comes before D/L
    #$name =~ s/(?<!\w)(d|l|dl)(?!\w)(.*?)(?<!\w)\((\+|\-|\+\-)\)(?!\w)/$1$2/g;
    #$name =~ s/(?<!\w)\((\+|\-|\+\-)\)(?!\w)(.*?)(?<!\w)(d|l|dl)(?!\w)/$3$2/g;
    $name =~ s/(?<!\w)(d|l|dl)-*\((\+|\-|\+\-)\)(?!\w)/$1/g;
    $name =~ s/(?<!\w)\((\+|\-|\+\-)\)-*(d|l|dl)(?!\w)/$2/g;

    ## protect -) as in (+/-) or (-) using capital letters    
    $name =~ s/[\(\[\{]-|-[\)\]\}]/MINUS/g;       # protect minus signs
    
    
    ## condense everything that isn't a letter, number, comma, or plus/minus
    ## except when between two numbers
    ##
    $name =~ s/[^,+\w]/-/g;           # convert to hyphens
    $name =~ s/-+/-/g;                # condense multiple hyphens in a row
    $name =~ s/(^-|-$)//g;            # strip leading/trailing hyphens
    $name =~ s/([,+])-/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/-([,+])/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/([a-z])-/$1/g;         # strip hyphens next to letters
    $name =~ s/-([a-z])/$1/g;         # strip hyphens near to letters
    $name =~ s/,+/,/g;                # condense multiple commas in a row
    $name =~ s/\++/\+/g;              # condense multiple pluses in a row
    
    ## de-protect and condense minus signs
    ##
    ## hypothetical example: 1-(-)-galacturonate   -->   1-galacturonate
    ##                       1-(+)-galacturonate   -->   1+galacturonate
    ##                       1-galacturonate       -->   1galacturonate
    $name =~ s/MINUS/-/g;
    $name =~ s/-+/-/g;

    # check for tandem duplicate names after conforming
    $name_len = length $name;
    if ($name_len % 2 == 0)
    {
        $half_name_len = 0.5 * $name_len;
        $name_half1 = substr $name, 0, $half_name_len;
        $name_half2 = substr $name, $half_name_len, $half_name_len;
        
        if ($name_half1 eq $name_half2)
        {
            $name = $name_half1;
        }
    }

    return $name;
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


sub cmp_merged_id
{
    my $id_1 = $a;
    my $id_2 = $b;
    my $formula_a  = '';
    my $formula_b  = '';
    my $name_a     = '';
    my $name_b     = '';
    my $inchi_a    = '';
    my $inchi_b    = '';
    my $hmdb_a     = '';
    my $hmdb_b     = '';
    my $kegg_a     = '';
    my $kegg_b     = '';
    my $cid_a      = '';
    my $cid_b      = '';
    my $lmaps_a    = '';
    my $lmaps_b    = '';
    my $compare    = '';
    my $name_len_a = '';
    my $name_len_b = '';
    
    my $name_a_conformed = '';
    my $name_b_conformed = '';


    # HMDB and PubChem have a lot of environmental metabolites,
    # put anything with a KEGG ID first, since we may care about it more
    if ($id_1 =~ /^HMDB/i)
    {
        $kegg_a = $hmdb_ccs_hash{$id_1}{kegg};
    }
    elsif ($id_1 =~ /^LM/i)
    {
        $kegg_a = $lmaps_ccs_hash{$id_1}{kegg};
    }
    if ($id_2 =~ /^HMDB/i)
    {
        $kegg_b = $hmdb_ccs_hash{$id_2}{kegg};
    }
    elsif ($id_2 =~ /^LM/i)
    {
        $kegg_b = $lmaps_ccs_hash{$id_2}{kegg};
    }
    if ($kegg_a ne '' && $kegg_b eq '') { return -1; }
    if ($kegg_a eq '' && $kegg_b ne '') { return  1; }
    

    # for now, don't sort on formula first
    if (0)
    {
        # first, sort on formula alphanumerically
        if ($id_1 =~ /^HMDB/i)
        {
            $formula_a = $hmdb_ccs_hash{$id_1}{formula};
        }
        elsif ($id_1 =~ /^LM/i)
        {
            $formula_a = $lmaps_ccs_hash{$id_1}{formula};
        }
        if ($id_2 =~ /^HMDB/i)
        {
            $formula_b = $hmdb_ccs_hash{$id_2}{formula};
        }
        elsif ($id_2 =~ /^LM/i)
        {
            $formula_b = $lmaps_ccs_hash{$id_2}{formula};
        }
        if ($formula_a ne '' && $formula_b eq '') { return -1; }
        if ($formula_a eq '' && $formula_b ne '') { return  1; }
        $compare = cmp_args_alphanumeric($formula_a, $formula_b);
        if ($compare) { return $compare; }
    }


    # then on name, alphanumerically
    if ($id_1 =~ /^HMDB/i)
    {
        $name_a = $hmdb_ccs_hash{$id_1}{name};
    }
    elsif ($id_1 =~ /^LM/i)
    {
        $name_a = $lmaps_ccs_hash{$id_1}{name};
    }
    if ($id_2 =~ /^HMDB/i)
    {
        $name_b = $hmdb_ccs_hash{$id_2}{name};
    }
    elsif ($id_2 =~ /^LM/i)
    {
        $name_b = $lmaps_ccs_hash{$id_2}{name};
    }
    if ($name_a ne '' && $name_b eq '') { return -1; }
    if ($name_a eq '' && $name_b ne '') { return  1; }
    
    $name_a_conformed = conform_name($name_a);
    $name_b_conformed = conform_name($name_b);


    # sort on conformed name length and name
    # then on original length and name
    # shorter names first
    #
    $name_len_a = length $name_a_conformed;
    $name_len_b = length $name_b_conformed;
    $compare = $name_len_a <=> $name_len_b;
    if ($compare) { return $compare; }

    $compare = cmp_args_alphanumeric($name_a_conformed, $name_b_conformed);
    if ($compare) { return $compare; }

    $name_len_a = length $name_a;
    $name_len_b = length $name_b;
    $compare = $name_len_a <=> $name_len_b;
    if ($compare) { return $compare; }

    $compare = cmp_args_alphanumeric($name_a, $name_b);
    if ($compare) { return $compare; }
    
    
    # for now, just sort on identifier, don't want to code the rest up yet
    # formula and name should be good enough for now
    return ($id_1 cmp $id_2);
}



# begin main()

$max_hits       = 15;       # metabolites with > max hits remain unannotated
$mz_tol_ppm_min =  5;       # 10 PPM
$mz_tol_ppm_max = 20;       # 20 PPM
$mz_tol_ppm_inc =  5;       #  5 PPM
$ccs_frac_tol   =  0.03;
#$ccs_delta_tol = 16;       # 16 Angstroms^2


$mass_proton = 1.0072764665789;       # mass of proton,   2022 CODATA
$mass_e      = 0.0005485799090441;    # mass of electron, 2022 CODATA

# all adducts listed here are +1 or -1 charge
# Fein Lab has outdated calculations from values circa ~2003
#
# newer atomic weights taken from:
#   https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=X
#   where X is the symbol of the element
#
$adduct_offset_hash{'[M+H]+'}  =  $mass_proton;
$adduct_offset_hash{'[M+Na]+'} =  22.9897692820 - $mass_e;  # (Fein: 22.989218)
$adduct_offset_hash{'[M+K]+'}  =  38.9637064864 - $mass_e;  # (Fein: 38.963158)
$adduct_offset_hash{'[M-H]-'}  = -$mass_proton;
$adduct_offset_hash{'[M+Cl]-'} =  34.968852682  + $mass_e;  # (Fein: 34.969402)
$adduct_charge_hash{'[M+H]+'}  =  +1;
$adduct_charge_hash{'[M+Na]+'} =  +1;
$adduct_charge_hash{'[M+K]+'}  =  +1;
$adduct_charge_hash{'[M-H]-'}  =  -1;
$adduct_charge_hash{'[M+Cl]-'} =  -1;


$more_hmdb_annotation_headers[0] = 'super_class';
$more_hmdb_annotation_headers[1] = 'class';
$more_hmdb_annotation_headers[2] = 'sub_class';
$more_hmdb_annotation_headers[3] = 'direct_parent';
$more_hmdb_annotation_headers[4] = 'molecular_framework';

$more_lmaps_annotation_headers[0] = 'ABBREVIATION';
$more_lmaps_annotation_headers[1] = 'CATEGORY';
$more_lmaps_annotation_headers[2] = 'MAIN_CLASS';
$more_lmaps_annotation_headers[3] = 'SUB_CLASS';
$more_lmaps_annotation_headers[4] = 'CLASS_LEVEL4';


$hmdb_ccs_filename  = shift;
$lmaps_ccs_filename = shift;
$scils_ccs_filename = shift;
$cardinal_filename  = shift;    # files derived from cardinal pipeline output

if (!defined($hmdb_ccs_filename) || !defined($lmaps_ccs_filename) ||
    !defined($scils_ccs_filename))
{
    $program_name = basename($0);

    print STDERR "Usage: $program_name hmdb_ccs.txt lmaps_ccs.txt scils_ccs.txt [hcdist_mz.txt]\n";
    print STDERR "  Input files:\n";
    print STDERR "    parsed HMDB annotation with CCS predictions\n";
    print STDERR "    parsed LipidMaps annotation with CCS predictions\n";
    print STDERR "    SCiLS exported m/z and CCS values for .imzML file\n";
    print STDERR "    m/z clusters from CARDINAL + hcdist pipeline\n";
    print STDERR "\n";
    print STDERR "  If the 4th file is not specified, output annotated SCiLS instead\n";
    
    exit(1);
}

$scils_file_present = 1;

# no SCiLS file provided, try to make do without it
#
# TODO -- attempt to improve guessing of ion mode some more
#
# Still guesses wrong ion mode too often, leave disabled for now
#
if (0 && $scils_ccs_filename =~ /^NULL$/i)
{
    print STDERR "WARNING -- no SCiLS file provided, guessing ion mode and adducts\n";

    $scils_file_present = 0;
}

if ($scils_file_present)
{
    open SCILS_CCS, "$scils_ccs_filename" or die "ABORT -- cannot open file $scils_ccs_filename\n";
}
open HMDB_CCS,  "$hmdb_ccs_filename"  or die "ABORT -- cannot open file $hmdb_ccs_filename\n";
open LMAPS_CCS, "$lmaps_ccs_filename" or die "ABORT -- cannot open file $lmaps_ccs_filename\n";

if (defined($cardinal_filename))
{
    open CARDINAL,  "$cardinal_filename"  or die "ABORT -- cannot open file $cardinal_filename\n";
}


# read in SCiLS ccs file
if ($scils_file_present)
{
  # jump down to header line
  while(defined($line=<SCILS_CCS>))
  {
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    # found likely header line
    if ($c =~ /\S/)
    {
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

    if (!defined($scils_ccs_mz_col) &&
        $field =~ /^m\/?z\b/i)
    {
        $scils_ccs_mz_col = $i;
    }
    elsif (!defined($scils_ccs_mz_width_col) &&
           $field =~ /Interval Width/i)
    {
        $scils_ccs_mz_width_col = $i;
    }
    elsif (!defined($scils_ccs_ccs_col) &&
           $field =~ /^CCS\s*\[/i)
    {
        $scils_ccs_ccs_col = $i;
    }
    elsif (!defined($scils_ccs_ccs_width_col) &&
           $field =~ /^CCS Interval Width/i)
    {
        $scils_ccs_ccs_width_col = $i;
    }
    elsif (!defined($scils_ccs_mass_neutral_col) &&
           $field =~ /^Neutral\s*Mass/i)
    {
        $scils_ccs_mass_neutral_col = $i;
    }
    elsif (!defined($scils_ccs_adduct_col) &&
           $field =~ /^Notation$/i)
    {
        $scils_ccs_adduct_col = $i;
    }
  }
  #$scils_ccs_header_line = join "\t", @array;


  # read in the SCiLS ccs data
  while(defined($line=<SCILS_CCS>))
  {
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    @array = split /\t/, $line, -1;
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
    }

    $mz           = $array[$scils_ccs_mz_col];
    $mz_width     = $array[$scils_ccs_mz_width_col];
    $ccs          = $array[$scils_ccs_ccs_col];
    $ccs_width    = $array[$scils_ccs_ccs_width_col];
    $mass_neutral = $array[$scils_ccs_mass_neutral_col];
    $adduct       = $array[$scils_ccs_adduct_col];

    if (is_number($mz))
    {    
        $line_new = join "\t", @array;
        
        $scils_ccs_hash{$mz}{mz_width}     = $mz_width;
        $scils_ccs_hash{$mz}{ccs}          = $ccs;
        $scils_ccs_hash{$mz}{ccs_width}    = $ccs_width;
        $scils_ccs_hash{$mz}{mass_neutral} = $mass_neutral;
        $scils_ccs_hash{$mz}{adduct}       = $adduct;
        $scils_ccs_hash{$mz}{line}         = $line_new;
    }
    
    #printf "%s  %s  %s  %s  %s  %s\n",
    #    $mz, $mz_width, $ccs, $ccs_width, $mass_neutral, $adduct;
  }
  close SCILS_CCS;
}


# read in HMDB ccs file
# jump down to header line
while(defined($line=<HMDB_CCS>))
{
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    # found likely header line
    if ($c =~ /\S/)
    {
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

    if (!defined($hmdb_ccs_id_col) &&
        $field =~ /^accession/i)
    {
        $hmdb_ccs_id_col = $i;
    }
    elsif (!defined($hmdb_ccs_mass_col) &&
           $field =~ /mono_mass/i)
    {
        $hmdb_ccs_mass_col = $i;
    }
    elsif (!defined($hmdb_ccs_charge_formal_col) &&
           $field =~ /formal_charge/i)
    {
        $hmdb_ccs_charge_formal_col = $i;
    }
    elsif (!defined($hmdb_ccs_charge_physio_col) &&
           $field =~ /physiological_charge/i)
    {
        $hmdb_ccs_charge_physio_col = $i;
    }
    elsif (!defined($hmdb_ccs_name_col) &&
           $field =~ /^name$/i)
    {
        $hmdb_ccs_name_col = $i;
    }
    elsif (!defined($hmdb_ccs_inchi_col) &&
           $field =~ /^inchi$/i)
    {
        $hmdb_ccs_inchi_col = $i;
    }
    elsif (!defined($hmdb_ccs_kegg_col) &&
           $field =~ /^kegg_id/i)
    {
        $hmdb_ccs_kegg_col = $i;
    }
    elsif (!defined($hmdb_ccs_cid_col) &&
           $field =~ /^pubchem_compound_id/i)
    {
        $hmdb_ccs_cid_col = $i;
    }
    
    # adduct CCS
    if ($field =~ /^CCS:\s*(.*)/)
    {
        $adduct = $1;

        $hmdb_adduct_col_hash{$adduct} = $i;
    }
    
    # formula columns
    if (!defined($hmdb_ccs_formula_orig_col) &&
        $field =~ /^chemical_formula/i)
    {
        $hmdb_ccs_formula_orig_col = $i;
    }
    elsif (!defined($hmdb_ccs_formula_chosen_col) &&
        $field =~ /^formula_chosen/i)
    {
        $hmdb_ccs_formula_chosen_col = $i;
    }
    
    $hmdb_header_col_hash{$field} = $i;
}
#$hmdb_ccs_header_line = join "\t", @array;


# prefer formula_chosen columns over original chemical_formula column
if (defined($hmdb_ccs_formula_chosen_col))
{
    $hmdb_ccs_formula_col = $hmdb_ccs_formula_chosen_col;
}
elsif (defined($hmdb_ccs_formula_orig_col))
{
    $hmdb_ccs_formula_col = $hmdb_ccs_formula_orig_col;
}


# scan for additional annotation columns
for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
{
    $header = $more_hmdb_annotation_headers[$i];
    $col    = $hmdb_header_col_hash{$header};

    if (defined($col))
    {
        $more_hmdb_annotation_cols[$i] = $col;
    }
}


# read in the HMDB ccs data
while(defined($line=<HMDB_CCS>))
{
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    @array = split /\t/, $line, -1;
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
    }

    $hmdb_id            = $array[$hmdb_ccs_id_col];
    $hmdb_mass          = $array[$hmdb_ccs_mass_col];
    $hmdb_formula       = $array[$hmdb_ccs_formula_col];
    $hmdb_charge_formal = $array[$hmdb_ccs_charge_formal_col];
    $hmdb_charge_physio = $array[$hmdb_ccs_charge_physio_col];
    $hmdb_name          = $array[$hmdb_ccs_name_col];
    $hmdb_inchi         = $array[$hmdb_ccs_inchi_col];
    $hmdb_kegg          = $array[$hmdb_ccs_kegg_col];
    $hmdb_cid           = $array[$hmdb_ccs_cid_col];
    
    $hmdb_ccs_mph       = $array[$hmdb_adduct_col_hash{'[M+H]+'}];
    $hmdb_ccs_mpna      = $array[$hmdb_adduct_col_hash{'[M+Na]+'}];
    $hmdb_ccs_mmh       = $array[$hmdb_adduct_col_hash{'[M-H]-'}];

    # skip unexpected errors
    if (!is_number($hmdb_mass)) { next; }

    
    $hmdb_ccs_hash{$hmdb_id}{mass}          = $hmdb_mass;
    $hmdb_ccs_hash{$hmdb_id}{formula}       = $hmdb_formula;
    $hmdb_ccs_hash{$hmdb_id}{charge_formal} = $hmdb_charge_formal;
    $hmdb_ccs_hash{$hmdb_id}{charge_physio} = $hmdb_charge_physio;
    $hmdb_ccs_hash{$hmdb_id}{name}          = $hmdb_name;
    $hmdb_ccs_hash{$hmdb_id}{inchi}         = $hmdb_inchi;
    $hmdb_ccs_hash{$hmdb_id}{kegg}          = $hmdb_kegg;
    $hmdb_ccs_hash{$hmdb_id}{cid}           = $hmdb_cid;

    if (is_number($hmdb_ccs_mph))
    {
        $hmdb_ccs_hash{$hmdb_id}{'[M+H]+'}      = $hmdb_ccs_mph;
        $hmdb_ccs_hash{$hmdb_id}{has_ccs}       = 1;
    }
    if (is_number($hmdb_ccs_mpna))
    {
        $hmdb_ccs_hash{$hmdb_id}{'[M+Na]+'}     = $hmdb_ccs_mpna;
        $hmdb_ccs_hash{$hmdb_id}{has_ccs}       = 1;
    }
    if (is_number($hmdb_ccs_mmh))
    {
        $hmdb_ccs_hash{$hmdb_id}{'[M-H]-'}      = $hmdb_ccs_mmh;
        $hmdb_ccs_hash{$hmdb_id}{has_ccs}       = 1;
    }


    # read in additional annotations
    for ($i = 0; $i < @more_hmdb_annotation_cols; $i++)
    {
        $col = $more_hmdb_annotation_cols[$i];
        
        if (defined($col))
        {
            $more_hmdb_annotation_hash{$hmdb_id}[$i] = $array[$col];
        }
    }

    
    #printf "%s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n",
    #    $hmdb_id, $hmdb_mass, $hmdb_formula,
    #    $hmdb_charge_formal, $hmdb_charge_physio,
    #    $hmdb_name, $hmdb_inchi, $hmdb_kegg, $hmdb_cid,
    #    $hmdb_ccs_mph, $hmdb_ccs_mpna, $hmdb_ccs_mmh;
    
}
close HMDB_CCS;



# read in LMAPS ccs file
# jump down to header line
while(defined($line=<LMAPS_CCS>))
{
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    # found likely header line
    if ($c =~ /\S/)
    {
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

    if (!defined($lmaps_ccs_id_col) &&
        $field =~ /^LM_ID/i)
    {
        $lmaps_ccs_id_col = $i;
    }
    elsif (!defined($lmaps_ccs_mass_col) &&
           $field =~ /^EXACT_MASS/i)
    {
        $lmaps_ccs_mass_col = $i;
    }
    elsif (!defined($lmaps_ccs_name_col) &&
           $field =~ /^CHOSEN_NAME/i)
    {
        $lmaps_ccs_name_col = $i;
    }
    elsif (!defined($lmaps_ccs_inchi_col) &&
           $field =~ /^INCHI$/i)
    {
        $lmaps_ccs_inchi_col = $i;
    }
    elsif (!defined($lmaps_ccs_cid_col) &&
           $field =~ /^PUBCHEM_CID/i)
    {
        $lmaps_ccs_cid_col = $i;
    }
    elsif (!defined($lmaps_ccs_hmdb_col) &&
           $field =~ /^HMDB_ID/i)
    {
        $lmaps_ccs_hmdb_col = $i;
    }
    elsif (!defined($lmaps_ccs_kegg_col) &&
           $field =~ /^KEGG_ID/i)
    {
        $lmaps_ccs_kegg_col = $i;
    }

    # adduct CCS
    if ($field =~ /^CCS:\s*(.*)/)
    {
        $adduct = $1;

        $lmaps_adduct_col_hash{$adduct} = $i;
    }
    
    # formula columns
    if (!defined($lmaps_ccs_formula_orig_col) &&
        $field =~ /^FORMULA/i)
    {
        $lmaps_ccs_formula_orig_col = $i;
    }
    elsif (!defined($lmaps_ccs_formula_chosen_col) &&
        $field =~ /^formula_chosen/i)
    {
        $lmaps_ccs_formula_chosen_col = $i;
    }
    
    $lmaps_header_col_hash{$field} = $i;
}
#$lmaps_ccs_header_line = join "\t", @array;


# prefer formula_chosen columns over original chemical_formula column
if (defined($lmaps_ccs_formula_chosen_col))
{
    $lmaps_ccs_formula_col = $lmaps_ccs_formula_chosen_col;
}
elsif (defined($lmaps_ccs_formula_orig_col))
{
    $lmaps_ccs_formula_col = $lmaps_ccs_formula_orig_col;
}


# scan for additional annotation columns
for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
{
    $header = $more_lmaps_annotation_headers[$i];
    $col    = $lmaps_header_col_hash{$header};

    if (defined($col))
    {
        $more_lmaps_annotation_cols[$i] = $col;
    }
}


# read in the LMAPS ccs data
while(defined($line=<LMAPS_CCS>))
{
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    @array = split /\t/, $line, -1;
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
    }

    $lmaps_id            = $array[$lmaps_ccs_id_col];
    $lmaps_mass          = $array[$lmaps_ccs_mass_col];
    $lmaps_formula       = $array[$lmaps_ccs_formula_col];
    $lmaps_name          = $array[$lmaps_ccs_name_col];
    $lmaps_inchi         = $array[$lmaps_ccs_inchi_col];
    $lmaps_cid           = $array[$lmaps_ccs_cid_col];
    $lmaps_hmdb          = $array[$lmaps_ccs_hmdb_col];
    $lmaps_kegg          = $array[$lmaps_ccs_kegg_col];
    
    $lmaps_ccs_mph       = $array[$lmaps_adduct_col_hash{'[M+H]+'}];
    $lmaps_ccs_mpna      = $array[$lmaps_adduct_col_hash{'[M+Na]+'}];
    $lmaps_ccs_mmh       = $array[$lmaps_adduct_col_hash{'[M-H]-'}];

    # skip unexpected errors
    if (!is_number($lmaps_mass)) { next; }

    # sanity check HMDB mapping (example: HMS3650A16)
    if (!($lmaps_hmdb =~ /^HMDB/i))
    {
        $lmaps_hmdb = '';
    }
    
    $lmaps_ccs_hash{$lmaps_id}{mass}    = $lmaps_mass;
    $lmaps_ccs_hash{$lmaps_id}{formula} = $lmaps_formula;
    $lmaps_ccs_hash{$lmaps_id}{name}    = $lmaps_name;
    $lmaps_ccs_hash{$lmaps_id}{inchi}   = $lmaps_inchi;
    $lmaps_ccs_hash{$lmaps_id}{cid}     = $lmaps_cid;
    $lmaps_ccs_hash{$lmaps_id}{hmdb}    = $lmaps_hmdb;
    $lmaps_ccs_hash{$lmaps_id}{kegg}    = $lmaps_kegg;

    if (is_number($lmaps_ccs_mph))
    {
        $lmaps_ccs_hash{$lmaps_id}{'[M+H]+'}      = $lmaps_ccs_mph;
        $lmaps_ccs_hash{$lmaps_id}{has_ccs}       = 1;
    }
    if (is_number($lmaps_ccs_mpna))
    {
        $lmaps_ccs_hash{$lmaps_id}{'[M+Na]+'}     = $lmaps_ccs_mpna;
        $lmaps_ccs_hash{$lmaps_id}{has_ccs}       = 1;
    }
    if (is_number($lmaps_ccs_mmh))
    {
        $lmaps_ccs_hash{$lmaps_id}{'[M-H]-'}      = $lmaps_ccs_mmh;
        $lmaps_ccs_hash{$lmaps_id}{has_ccs}       = 1;
    }


    # read in additional annotations
    for ($i = 0; $i < @more_lmaps_annotation_cols; $i++)
    {
        $col = $more_lmaps_annotation_cols[$i];
        
        if (defined($col))
        {
            $more_lmaps_annotation_hash{$lmaps_id}[$i] = $array[$col];
        }
    }


    #printf "%s  %s  %s  %s  %s  %s  %s  %s  %s %s  %s\n",
    #    $lmaps_id, $lmaps_mass, $lmaps_formula, $lmaps_name,
    #    $lmaps_inchi, $lmaps_cid, $lmaps_hmdb, $lmaps_kegg,
    #    $lmaps_ccs_mph, $lmaps_ccs_mpna, $lmaps_ccs_mmh
}
close LMAPS_CCS;


# read in CARDINAL ccs file
if (defined($cardinal_filename))
{
  # jump down to header line
  while(defined($line=<CARDINAL>))
  {
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    # found likely header line
    if ($c =~ /\S/)
    {
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

    if (!defined($cardinal_mz_col) &&
        $field =~ /^m\/?z/i)
    {
        $cardinal_mz_col = $i;
    }
  }
  if (!defined($cardinal_mz_col))
  {
    $cardinal_mz_col = 0;
  }
  $cardinal_header_line = join "\t", @array;


  # read in the CARDINAL ccs data
  $row = 0;
  while(defined($line=<CARDINAL>))
  {
    # skip comment lines
    $c = substr $line, 0, 1;
    if ($c eq '#') { next; }

    @array = split /\t/, $line, -1;
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
    }

    $cardinal_mz = $array[$cardinal_mz_col];
    $cardinal_mz =~ s/^[^_]+_//;
    
    if (is_number($cardinal_mz))
    {
        $line_new                            = join "\t", @array;
        $cardinal_row_array[$row]{line}      = $line_new;
        $cardinal_row_array[$row]{mz}        = $cardinal_mz;
        $cardinal_seen_mz_hash{$cardinal_mz} = 1;

        $row++;
    }
  }
  close CARDINAL;
  
  #$num_rows_cardinal = $row;
}

@merged_id_array        = sort (keys %hmdb_ccs_hash, keys %lmaps_ccs_hash);
@supported_adduct_array = sort keys %adduct_offset_hash;

foreach $merged_id (@merged_id_array)
{
    $mass          = '';
    $formal_charge = 0;

    if ($merged_id =~ /^HMDB/)
    {
        $mass          = $hmdb_ccs_hash{$merged_id}{mass};
        $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

        if (!defined($formal_charge) || !is_number($formal_charge))
        {
            $formal_charge = 0;
        }
    }
    elsif ($merged_id =~ /^LM/)
    {
        $mass = $lmaps_ccs_hash{$merged_id}{mass};
    }
    
    # skip on unexpected errors
    if (!is_number($mass)) { next; }
    
    # skip formal charges, since they are more likely in LC than MALDI
    if ($formal_charge)    { next; }

    
    # bin m/z for rapid scanning later
    # +/-1 is ~1000 ppm, which is *way* over our 10 ppm tolerance
    # bins are too large, but are very simple to code up
    #
    # also factor in +/- 2.014552 (pos - neg) difference
    #
    foreach $adduct (@supported_adduct_array)
    {
        $adduct_offset = $adduct_offset_hash{$adduct};
        
        # deal with formal charges
        if ($formal_charge)
        {
            $adduct_charge = $adduct_charge_hash{$adduct};

            # can never switch between opposite charges, skip it
            if (($formal_charge > 0 && $adduct_charge < 0) ||
                ($formal_charge < 0 && $adduct_charge > 0))
            {
                next;
            }
        
            # zero out H adduct offsets
            if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
            {
                $adduct_offset = 0;
            }
        
            $delta_charge   = $adduct_charge - $formal_charge;
            $adduct_offset += $delta_charge * $mass_proton;
        }

        $mz_new     = $mass + $adduct_offset;

        $ppm_offset = $mz_tol_ppm_max * $mz_new / 1000000;

        $mz_floor   = floor($mz_new);
        $mz_minus1  = floor($mz_new - $ppm_offset);
        $mz_plus1   = floor($mz_new + $ppm_offset);
    
        $mz_id_bins_hash{$mz_minus1}{$merged_id} = 1;
        $mz_id_bins_hash{$mz_floor}{$merged_id}  = 1;
        $mz_id_bins_hash{$mz_plus1}{$merged_id}  = 1;
    }
}


@scils_mz_array = sort {$a <=> $b} keys %scils_ccs_hash;

$i = 0;
for ($ppm = $mz_tol_ppm_min; $ppm <= $mz_tol_ppm_max; $ppm += $mz_tol_ppm_inc)
{
    $ppm_tol_array[$i++] = $ppm;
}


@mz_matched_array  = sort {$a <=> $b} keys %scils_match_hash;
$count_matched_neg = 0;
$count_matched_pos = 0;


# determine which ion mode we ran in, based on annotated adducts
if ($scils_file_present)
{
  foreach $mz_scils (@scils_mz_array)
  {
    $adduct = $scils_ccs_hash{$mz_scils}{adduct};
    
    # skip unsupported adducts
    if (!defined($adduct_offset_hash{$adduct})) { next; }
    
    if ($adduct =~ /\]\-/)
    {
        $count_matched_neg++;
    }
    if ($adduct =~ /\]\+/)
    {
        $count_matched_pos++;
    }
  }
  printf STDERR "#POS\t%s\t#NEG\t%s\n", $count_matched_pos, $count_matched_neg;
}
# uh oh, try to guess the ion mode...
else
{
    %trial_adduct_match_hash = ();
    %trial_mz_match_hash     = ();

    $trial_adduct_array[0] = '[M+H]+';
    $trial_adduct_array[1] = '[M-H]-';
    #@trial_adduct_array = @supported_adduct_array;

    # strict ppm tol for ion mode guesses
    #@trial_ppm_tol_array = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
    @trial_ppm_tol_array = ( 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
    #@trial_ppm_tol_array = ( 5, 10, 15, 20 );


    # m/z from cardinal reprocessing, etc.
    @cardinal_mz_array = sort {$a <=> $b} keys %cardinal_seen_mz_hash;
    
    foreach $cardinal_mz (@cardinal_mz_array)
    {
        $mz_floor = floor($cardinal_mz);

        foreach $ppm_tol (@trial_ppm_tol_array)
        {
            $found_match_flag = 0;
            foreach $adduct (@trial_adduct_array)
            {
                $count = 0;

                # will include MANY outside the acceptable tolerance,
                # but is quick and dirty to code up
                #
                @crude_candidate_array = ();
                if (defined($mz_id_bins_hash{$mz_floor}))
                {
                    @crude_candidate_array =
                        sort keys %{$mz_id_bins_hash{$mz_floor}};
                }
                
                # refine crude candidate list by m/z tolerance
                foreach $merged_id (@crude_candidate_array)
                {
                    $mass_db       = '';
                    $formal_charge = 0;
                    
                    if ($merged_id =~ /^HMDB/)
                    {
                        $mass_db       = $hmdb_ccs_hash{$merged_id}{mass};
                        $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

                        if (!defined($formal_charge) || !is_number($formal_charge))
                        {
                            $formal_charge = 0;
                        }
                    }
                    elsif ($merged_id =~ /^LM/)
                    {
                        $mass_db = $lmaps_ccs_hash{$merged_id}{mass};
                        
                        # skip LipidMaps during ion mode guessing
                        # don't want to deal with de-duplication
                        # HMDB has a lot of lipids anyways
                        next;
                    }
                    

                    if (!is_number($mass_db)) { next; }

                    # skip formal charges, since they are more likely in LC than MALDI
                    if ($formal_charge) { next; }

                    $adduct_offset = $adduct_offset_hash{$adduct};
                    if (!defined($adduct_offset)) { next; }

                    # deal with formal charges
                    if ($formal_charge)
                    {
                        $adduct_charge = $adduct_charge_hash{$adduct};

                        # can never switch between opposite charges, skip it
                        if (($formal_charge > 0 && $adduct_charge < 0) ||
                            ($formal_charge < 0 && $adduct_charge > 0))
                        {
                            next;
                        }
                    
                        # zero out H adduct offsets
                        if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
                        {
                            $adduct_offset = 0;
                        }
                    
                        $delta_charge   = $adduct_charge - $formal_charge;
                        $adduct_offset += $delta_charge * $mass_proton;
                    }

                    
                    $mz_db = $mass_db + $adduct_offset;
                    $ppm   = 1000000 * abs(($mz_db - $cardinal_mz) / $mz_db);
                    
                    if ($ppm <= $ppm_tol + 1.0e-5)
                    {
                        $trial_adduct_match_hash{$adduct}{$cardinal_mz}{merged_id}{$merged_id}
                            = 1;
                        $trial_adduct_match_hash{$adduct}{$cardinal_mz}{ppm_tol}
                            = $ppm_tol;
                        
                        $trial_mz_match_hash{$cardinal_mz}{$adduct} = 1;

                        #printf STDERR "TRIAL_FOOBAR\t%s\t%s\n",
                        #    $cardinal_mz, $merged_id;

                        $found_match_flag = 1;
                        
                        $count++;
                    }
                }

                if ($count)
                {
                    $trial_adduct_match_hash{$adduct}{$cardinal_mz}{count} = $count;
                }
            }

            if ($found_match_flag)
            {
                $ppm_tol_kept = $ppm_tol;
                last;
            }
        }
        
        # if mz matched multiple adducts, discard all the matches
        if (defined($trial_mz_match_hash{$cardinal_mz}))
        {
            @matched_adduct_array =
                keys %{$trial_mz_match_hash{$cardinal_mz}};
            
            if (@matched_adduct_array >= 2)
            {
                foreach $adduct (@matched_adduct_array)
                {
                    delete $trial_adduct_match_hash{$adduct}{$cardinal_mz};
                    delete $trial_mz_match_hash{$cardinal_mz}{$adduct};
                }
            }
            else
            {
                $adduct = $matched_adduct_array[0];
            
                $count =
                    $trial_adduct_match_hash{$adduct}{$cardinal_mz}{count};
                $ppm_tol =
                    $trial_adduct_match_hash{$adduct}{$cardinal_mz}{ppm_tol};

                # only keep somewhat unique hits
                if (1 || $count <= 15)
                {
                    #printf STDERR "TRIAL\t%s\t%.4f\t%s\t%s\n",
                    #    $adduct, $cardinal_mz, $ppm_tol, $count;
                }
                else
                {
                    delete $trial_adduct_match_hash{$adduct}{$cardinal_mz};
                    delete $trial_mz_match_hash{$cardinal_mz}{$adduct};
                }
            }
        }
    }
    
    foreach $adduct (@trial_adduct_array)
    {
        $sum         = 0.0;
        $avg_ppm_tol = 0.0;
    
        @matched_mz_array = ();
        
        if (defined($trial_adduct_match_hash{$adduct}))
        {
            @matched_mz_array = keys %{$trial_adduct_match_hash{$adduct}};
        }
        
        foreach $matched_mz (@matched_mz_array)
        {
            $count   = $trial_adduct_match_hash{$adduct}{$matched_mz}{count};
            $ppm_tol = $trial_adduct_match_hash{$adduct}{$matched_mz}{ppm_tol};
            
            $score  = 1.0 / ($ppm_tol);
            $score /= sqrt($count);
            #$score   = 1.0;
            
            $avg_ppm_tol += $ppm_tol;
            
            $sum += $score;
        }
        
        if (@matched_mz_array)
        {
            $avg_ppm_tol /= @matched_mz_array;
        }
        # don't divide by zero
        else
        {
            $avg_ppm_tol = 1.0;
        }
        
        $sum /= ($avg_ppm_tol * $avg_ppm_tol);

        if ($adduct =~ /\]\-/)
        {
            $count_matched_neg += $sum;
        }
        if ($adduct =~ /\]\+/)
        {
            $count_matched_pos += $sum;
        }

        
        printf STDERR "TRIAL_SCORE\t%s\t%s\t%s\n",
            $adduct, $avg_ppm_tol, $sum;
    }

    printf STDERR "#POS\t%s\t#NEG\t%s\n",
        $count_matched_pos, $count_matched_neg;
}



@assume_adduct_array = ();
$num_assumed_adduct_types = 0;
if ($count_matched_pos > $count_matched_neg)
{
    $assume_adduct_array[$num_assumed_adduct_types++] = '[M+H]+';
}
if ($count_matched_neg > $count_matched_pos)
{
    $assume_adduct_array[$num_assumed_adduct_types++] = '[M-H]-';
}


foreach $mz_scils (@scils_mz_array)
{
    $adduct = $scils_ccs_hash{$mz_scils}{adduct};

    # guess missing adduct
    if (!($adduct =~ /[A-Za-z]/) && $num_assumed_adduct_types)
    {
        $adduct = $assume_adduct_array[0];
        
        printf STDERR "Assuming missing SCiLS adduct:\t%s\t%s\n",
            $mz_scils, $adduct;
    }
    
    # skip unsupported adducts
    if (!defined($adduct_offset_hash{$adduct})) { next; }
    

    $mz_floor = floor($mz_scils);
    
    # crude bin lookup of potential candidate molecules;
    # will include MANY outside the acceptable tolerance,
    # but is quick and dirty to code up
    #
    @crude_candidate_array = ();
    if (defined($mz_id_bins_hash{$mz_floor}))
    {
        @crude_candidate_array = sort keys %{$mz_id_bins_hash{$mz_floor}};
    }
    
    #$count = @crude_candidate_array;
    #printf STDERR "%s\t%s\n", $mz_scils, $count;
    
    # refine crude candidate list by m/z tolerance
    %mz_candidate_lax_hash = ();
    foreach $merged_id (@crude_candidate_array)
    {
        $mass_db       = '';
        $formal_charge = 0;

        if ($merged_id =~ /^HMDB/)
        {
            $mass_db       = $hmdb_ccs_hash{$merged_id}{mass};
            $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

            if (!defined($formal_charge) || !is_number($formal_charge))
            {
                $formal_charge = 0;
            }
        }
        elsif ($merged_id =~ /^LM/)
        {
            $mass_db = $lmaps_ccs_hash{$merged_id}{mass};
        }
        

        if (!is_number($mass_db)) { next; }

        # skip formal charges, since they are more likely in LC than MALDI
        if ($formal_charge)       { next; }
        

        $adduct_offset = $adduct_offset_hash{$adduct};
        if (!defined($adduct_offset)) { next; }

        # deal with formal charges
        if ($formal_charge)
        {
            $adduct_charge = $adduct_charge_hash{$adduct};

            # can never switch between opposite charges, skip it
            if (($formal_charge > 0 && $adduct_charge < 0) ||
                ($formal_charge < 0 && $adduct_charge > 0))
            {
                next;
            }
        
            # zero out H adduct offsets
            if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
            {
                $adduct_offset = 0;
            }
        
            $delta_charge   = $adduct_charge - $formal_charge;
            $adduct_offset += $delta_charge * $mass_proton;
        }

        
        $mz_db = $mass_db + $adduct_offset;
        $ppm   = 1000000 * abs(($mz_db - $mz_scils) / $mz_db);
        
        if ($ppm <= $mz_tol_ppm_max + 1.0e-5)
        {
            $mz_candidate_lax_hash{$merged_id} = $ppm;
        }
    }
    @mz_candidate_lax_array = sort keys %mz_candidate_lax_hash;

    $count = @mz_candidate_lax_array;

    # no m/z hits within maximum ppm tolerance
    if ($count == 0) { next; }

    #printf STDERR "%s\t%s\n", $mz_scils, $count;


    # check for CCS matches at increasingly large ppm tolerances
    @kept_match_array         = ();
    $ppm_tol_kept             = 0;
    $matched_on_ccs_pass_flag = 0;
    $found_match_flag         = 0;

    foreach $ppm_tol (@ppm_tol_array)
    {
        %match_hash       = ();

        foreach $merged_id (@mz_candidate_lax_array)
        {
            $ppm = $mz_candidate_lax_hash{$merged_id};

            # skip ppm that are too large
            if ($ppm > $ppm_tol + 1e-5) { next; }

            $ccs_scils = $scils_ccs_hash{$mz_scils}{ccs};

            # invalid or missing CCS, skip it
            if (!defined($ccs_scils) || !is_number($ccs_scils))
            {
                next;
            }

            $ccs_db = '';
            if ($merged_id =~ /^HMDB/)
            {
                $ccs_db = $hmdb_ccs_hash{$merged_id}{$adduct};
            }
            elsif ($merged_id =~ /^LM/)
            {
                $ccs_db = $lmaps_ccs_hash{$merged_id}{$adduct};
            }

            $ccs_delta = 999;
            $ccs_frac  = 999;
            if (defined($ccs_db) && is_number($ccs_db))
            {
                $ccs_delta = abs($ccs_db - $ccs_scils);
                $ccs_frac  = $ccs_delta / $ccs_db;
            }

            #if ($ccs_delta <= $ccs_delta_tol + 1e-5 ||
            #    $ccs_frac  <= $ccs_frac_tol  + 1e-5)
            #if (defined($ccs_db) && is_number($ccs_db))

            if ($ccs_frac  <= $ccs_frac_tol  + 1e-5)
            {
                $match_hash{$merged_id} = $ppm;
                $found_match_flag       = 1;

                #if ($mz_scils > 168 && $mz_scils < 168.006)
                #{
                #    printf "%s  %s  %s  %s  %s  %s\n",
                #        $mz_scils, $merged_id,
                #        $ccs_scils, $ccs_db,
                #        $ccs_delta, $ccs_frac;
                #}
            }
        }
        
        if ($found_match_flag)
        {
            @kept_match_array         = sort keys %match_hash;
            $ppm_tol_kept             = $ppm_tol;
            $matched_on_ccs_pass_flag = 1;

            last;
        }
    }
    
    # 2nd pass, no CCS matches, only m/z
    # require no CCS predictions to exist
    if (@kept_match_array == 0)
    {
        foreach $ppm_tol (@ppm_tol_array)
        {
            %match_hash       = ();
            $found_match_flag = 0;

            foreach $merged_id (@mz_candidate_lax_array)
            {
                $ppm = $mz_candidate_lax_hash{$merged_id};

                # skip ppm that are too large
                if ($ppm > $ppm_tol + 1e-5) { next; }

                $has_ccs_flag = 0;
                if ($merged_id =~ /^HMDB/)
                {
                    if (defined($hmdb_ccs_hash{$merged_id}{has_ccs}))
                    {
                        $has_ccs_flag = 1;
                    }
                }
                elsif ($merged_id =~ /^LM/)
                {
                    if (defined($lmaps_ccs_hash{$merged_id}{has_ccs}))
                    {
                        $has_ccs_flag = 1;
                    }
                }

                # only match entries with no CCS predictions at all;
                # [M-H]- may be missing due to no hydrogens present
                # this is for metabolites that are too big and/or have
                # too complicated stereochemistry for DarkChem to predict
                #
                if ($has_ccs_flag == 0)
                {
                    $match_hash{$merged_id} = $ppm;
                    $found_match_flag       = 1;
                }
            }
            
            if ($found_match_flag)
            {
                @kept_match_array = sort keys %match_hash;
                $ppm_tol_kept = $ppm_tol;

                last;
            }
        }
    }
    
    
    # populate match hash from kept array
    %match_hash = ();
    if (@kept_match_array)
    {
        foreach $merged_id (@kept_match_array)
        {
            $match_hash{$merged_id} = 1;
        }
    }
    
    
    # keep only HMDB identifiers when LipidMaps is the same molecule
    # dedupe on inchi, hmdb, kegg, cid
    if (@kept_match_array)
    {
        %temp_lmaps_dupe_hash  = ();    # is LipidMaps covered by HMDB hit?
        %temp_hmdb_hash        = ();
        %temp_hmdb_inchi_hash  = ();
        %temp_hmdb_kegg_hash   = ();
        %temp_hmdb_cid_hash    = ();
        %temp_lmaps_hash       = ();
        %temp_lmaps_inchi_hash = ();
        %temp_lmaps_kegg_hash  = ();
        %temp_lmaps_cid_hash   = ();
        
        foreach $merged_id (@kept_match_array)
        {
            if ($merged_id =~ /^HMDB/)
            {
                $temp_hmdb_hash{$merged_id} = 1;
                
                $inchi = $hmdb_ccs_hash{$merged_id}{inchi};
                $kegg  = $hmdb_ccs_hash{$merged_id}{kegg};
                $cid   = $hmdb_ccs_hash{$merged_id}{cid};
                
                if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                {
                    $temp_hmdb_inchi_hash{$inchi} = 1;
                }
                if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                {
                    $temp_hmdb_kegg_hash{$kegg} = 1;
                }
                if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                {
                    $temp_hmdb_cid_hash{$cid} = 1;
                }
            }
            elsif ($merged_id =~ /^LM/)
            {
                $temp_lmaps_hash{$merged_id} = 1;

                $hmdb  = $lmaps_ccs_hash{$merged_id}{hmdb};
                $inchi = $lmaps_ccs_hash{$merged_id}{inchi};
                $kegg  = $lmaps_ccs_hash{$merged_id}{kegg};
                $cid   = $lmaps_ccs_hash{$merged_id}{cid};
                
                if (defined($hmdb) && $hmdb =~ /[A-Za-z0-9]/)
                {
                    $temp_lmaps_hmdb_hash{$merged_id}{$hmdb} = 1;
                }
                if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                {
                    $temp_lmaps_inchi_hash{$merged_id}{$inchi} = 1;
                }
                if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                {
                    $temp_lmaps_kegg_hash{$merged_id}{$kegg} = 1;
                }
                if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                {
                    $temp_lmaps_cid_hash{$merged_id}{$cid} = 1;
                }
            }
        }
        
        @temp_hmdb_array       = sort keys %temp_hmdb_hash;
        @temp_lmaps_array      = sort keys %temp_lmaps_hash;
        @temp_lmaps_dupe_array = ();

        foreach $lmaps_id (@temp_lmaps_array)
        {
            $original_lmaps_matches{$mz_scils}{$lmaps_id} =
                $matched_on_ccs_pass_flag;
        }

        if (@temp_hmdb_array && @temp_lmaps_array)
        {
            foreach $lmaps_id (@temp_lmaps_array)
            {
                @temp_array = sort keys %{$temp_lmaps_hmdb_hash{$lmaps_id}};
                foreach $identifier (@temp_array)
                {
                    if (defined($temp_hmdb_hash{$identifier}))
                    {
                        $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                        next;
                    }
                }

                @temp_array = sort keys %{$temp_lmaps_inchi_hash{$lmaps_id}};
                foreach $identifier (@temp_array)
                {
                    if (defined($temp_hmdb_inchi_hash{$identifier}))
                    {
                        $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                        next;
                    }
                }

                @temp_array = sort keys %{$temp_lmaps_kegg_hash{$lmaps_id}};
                foreach $identifier (@temp_array)
                {
                    if (defined($temp_hmdb_kegg_hash{$identifier}))
                    {
                        $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                        next;
                    }
                }

                @temp_array = sort keys %{$temp_lmaps_cid_hash{$lmaps_id}};
                foreach $identifier (@temp_array)
                {
                    if (defined($temp_hmdb_cid_hash{$identifier}))
                    {
                        $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                        next;
                    }
                }
            }
        }
        
        @temp_lmaps_dupe_array = sort keys %temp_lmaps_dupe_hash;
        foreach $lmaps_id (@temp_lmaps_dupe_array)
        {
            delete $match_hash{$lmaps_id};
        }
    }
    
    # re-build kept match array after deleting LipidMaps HMDB dupes
    @kept_match_array = sort keys %match_hash;

    
    $count = @kept_match_array;
    
    if ($count > 0)
    {
        foreach $merged_id (@kept_match_array)
        {
            $scils_match_hash{$mz_scils}{$merged_id}{ppm_tol} = $ppm_tol_kept;
            $scils_match_hash{$mz_scils}{$merged_id}{matched_on_ccs} =
                $matched_on_ccs_pass_flag;
        }

        #$match_str = join ";", @kept_match_array;
        
        #printf STDERR "%.4f\t%s\t%s\t%s\n",
        #    $mz_scils, $ppm_tol_kept, $matched_on_ccs_pass_flag,
        #    $count;
    }
}


# output vendor m/z table
if (!defined($cardinal_filename))
{
  printf "%s",   'm/z';
  printf "\t%s", 'Adduct';
  printf "\t%s", 'CCS';
  printf "\t%s", 'PPM tolerence';
  printf "\t%s", 'CCS matched';
  printf "\t%s", 'Number of matches';
  printf "\t%s", 'Formula';
  printf "\t%s", 'Name';
  printf "\t%s", 'InChI';
  printf "\t%s", 'HMDB';
  printf "\t%s", 'KEGG';
  printf "\t%s", 'PubChemCID';
  printf "\t%s", 'LipidMaps';    # will include HMDB dupes that were removed

  # insert additional annotation
  for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
  {
      $col = $more_hmdb_annotation_cols[$i];
      
      if (defined($col))
      {
          printf "\t%s", $more_hmdb_annotation_headers[$i];
      }
  }
  for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
  {
      $col = $more_lmaps_annotation_cols[$i];
      
      if (defined($col))
      {
          printf "\t%s", $more_lmaps_annotation_headers[$i];
      }
  }

  printf "\n";

  @mz_matched_array = sort {$a <=> $b} keys %scils_match_hash;

  foreach $mz_scils (@mz_matched_array)
  {
    %seen_formula_hash        = ();
    %seen_name_conformed_hash = ();
    %seen_inchi_hash          = ();
    %seen_hmdb_hash           = ();
    %seen_kegg_hash           = ();
    %seen_cid_hash            = ();

    %seen_more_hmdb_hash      = ();
    %seen_more_lmaps_hash     = ();

    $formula_str = '';
    $name_str    = '';
    $inchi_str   = '';
    $hmdb_str    = '';
    $kegg_str    = '';
    $cid_str     = '';
    $lmaps_str   = '';

    # initialize more annotation strings
    for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
    {
        $more_hmdb_str_array[$i]  = '';
    }
    for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
    {
        $more_lmaps_str_array[$i] = '';
    }
    
    $ppm_tol        = 0;
    $matched_on_ccs = 0;
    
    $ccs         = $scils_ccs_hash{$mz_scils}{ccs};
    $adduct      = $scils_ccs_hash{$mz_scils}{adduct};

    @merged_id_array = sort cmp_merged_id
        keys %{$scils_match_hash{$mz_scils}};
    
    $num_matches = @merged_id_array;
    
    for ($i = 0; $i < $num_matches && $i < $max_hits; $i++)
    {
        $merged_id = $merged_id_array[$i];
    
        if ($scils_match_hash{$mz_scils}{$merged_id}{ppm_tol})
        {
            $ppm_tol = $scils_match_hash{$mz_scils}{$merged_id}{ppm_tol};
        }
        if ($scils_match_hash{$mz_scils}{$merged_id}{matched_on_ccs})
        {
            $matched_on_ccs = 1;
        }
    
        $formula = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $formula = $hmdb_ccs_hash{$merged_id}{formula};
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $formula = $lmaps_ccs_hash{$merged_id}{formula};
        }
        if ($formula =~ /[A-Za-z0-9]/ &&
            !defined($seen_formula_hash{$formula}))
        {
            if ($formula_str eq '')
            {
                $formula_str = $formula;
            }
            else
            {
                $formula_str .= ' | ' . $formula;
            }
            $seen_formula_hash{$formula} = 1;
        }

        $name = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $name = $hmdb_ccs_hash{$merged_id}{name};
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $name = $lmaps_ccs_hash{$merged_id}{name};
        }
        $name_conformed = conform_name($name);
        if ($name =~ /[A-Za-z0-9]/ &&
            !defined($seen_name_conformed_hash{$name_conformed}))
        {
            if ($name_str eq '')
            {
                $name_str = $name;
            }
            else
            {
                $name_str .= ' | ' . $name;
            }
            $seen_name_conformed_hash{$name_conformed} = 1;
        }

        $inchi = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $inchi = $hmdb_ccs_hash{$merged_id}{inchi};
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $inchi = $lmaps_ccs_hash{$merged_id}{inchi};
        }
        if ($inchi =~ /[A-Za-z0-9]/ &&
            !defined($seen_inchi_hash{$inchi}))
        {
            if ($inchi_str eq '')
            {
                $inchi_str = $inchi;
            }
            else
            {
                $inchi_str .= ' | ' . $inchi;
            }
            $seen_inchi_hash{$inchi} = 1;
        }

        $hmdb = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $hmdb = $merged_id;
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $hmdb = $lmaps_ccs_hash{$merged_id}{hmdb};
        }
        if ($hmdb =~ /[A-Za-z0-9]/ &&
            !defined($seen_hmdb_hash{$hmdb}))
        {
            if ($hmdb_str eq '')
            {
                $hmdb_str = $hmdb;
            }
            else
            {
                $hmdb_str .= ' | ' . $hmdb;
            }
            $seen_hmdb_hash{$hmdb} = 1;
        }

        $kegg = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $kegg = $hmdb_ccs_hash{$merged_id}{kegg};
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $kegg = $lmaps_ccs_hash{$merged_id}{kegg};
        }
        if ($kegg =~ /[A-Za-z0-9]/ &&
            !defined($seen_kegg_hash{$kegg}))
        {
            if ($kegg_str eq '')
            {
                $kegg_str = $kegg;
            }
            else
            {
                $kegg_str .= ' | ' . $kegg;
            }
            $seen_kegg_hash{$kegg} = 1;
        }

        $cid = '';
        if ($merged_id =~ /^HMDB/i)
        {
            $cid = $hmdb_ccs_hash{$merged_id}{cid};
        }
        elsif ($merged_id =~ /^LM/i)
        {
            $cid = $lmaps_ccs_hash{$merged_id}{cid};
        }
        if ($cid =~ /[A-Za-z0-9]/ &&
            !defined($seen_cid_hash{$cid}))
        {
            if ($cid_str eq '')
            {
                $cid_str = $cid;
            }
            else
            {
                $cid_str .= ' | ' . $cid;
            }
            $seen_cid_hash{$cid} = 1;
        }
    }

    @lmaps_id_array = ();
    if (defined($original_lmaps_matches{$mz_scils}))
    {
        @lmaps_id_array = sort cmp_merged_id keys %{$original_lmaps_matches{$mz_scils}};
        $lmaps_str      = join ' | ', @lmaps_id_array;
    }

    # more annotations
    foreach $merged_id (@merged_id_array)
    {
        if (!defined($more_hmdb_annotation_hash{$merged_id}))
        {
            next;
        }

        for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
        {
            $value = $more_hmdb_annotation_hash{$merged_id}[$i];

            if (defined($value) && $value =~ /[A-Za-z0-9]/ &&
                !(defined($seen_more_hmdb_hash{$i}) &&
                  defined($seen_more_hmdb_hash{$i}{$value})))
            {
                if ($more_hmdb_str_array[$i] eq '')
                {
                    $more_hmdb_str_array[$i] = $value;
                }
                else
                {
                    $more_hmdb_str_array[$i] .= ' | ' . $value;
                }

                $seen_more_hmdb_hash{$i}{$value} = 1;
            }
        }
    }

    foreach $merged_id (@lmaps_id_array)
    {
        if (!defined($more_lmaps_annotation_hash{$merged_id}))
        {
            next;
        }

        for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
        {
            $value = $more_lmaps_annotation_hash{$merged_id}[$i];

            if (defined($value) && $value =~ /[A-Za-z0-9]/ &&
                !(defined($seen_more_lmaps_hash{$i}) &&
                  defined($seen_more_lmaps_hash{$i}{$value})))
            {
                if ($more_lmaps_str_array[$i] eq '')
                {
                    $more_lmaps_str_array[$i] = $value;
                }
                else
                {
                    $more_lmaps_str_array[$i] .= ' | ' . $value;
                }

                $seen_more_lmaps_hash{$i}{$value} = 1;
            }
        }
    }


    printf "%s",   $mz_scils;
    printf "\t%s", $adduct;
    printf "\t%s", $ccs;
    printf "\t%s", $ppm_tol;
    printf "\t%s", $matched_on_ccs;
    printf "\t%s", $num_matches;
    printf "\t%s", $formula_str;
    printf "\t%s", unicode_to_ascii($name_str);
    printf "\t%s", $inchi_str;
    printf "\t%s", $hmdb_str;
    printf "\t%s", $kegg_str;
    printf "\t%s", $cid_str;
    printf "\t%s", $lmaps_str;

    # insert additional annotation
    for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
    {
        $col = $more_hmdb_annotation_cols[$i];
        
        if (defined($col))
        {
            printf "\t%s", $more_hmdb_str_array[$i];
        }
    }
    for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
    {
        $col = $more_lmaps_annotation_cols[$i];
        
        if (defined($col))
        {
            printf "\t%s", $more_lmaps_str_array[$i];
        }
    }

    printf "\n";
  }
}


# re-match to cardinal m/z
if (defined($cardinal_filename))
{
    # vendor m/z matched to annotation
    @scils_mz_all_array     = sort {$a <=> $b} keys %scils_ccs_hash;
    @scils_mz_matched_array = sort {$a <=> $b} keys %scils_match_hash;

    # m/z from cardinal reprocessing, etc.
    @cardinal_mz_array = sort {$a <=> $b} keys %cardinal_seen_mz_hash;
    
    $scils_mz_matched_start = 0;
    $scils_mz_all_start     = 0;
    foreach $cardinal_mz (@cardinal_mz_array)
    {
        @kept_match_array          = ();
        %crude_mz_match_hash       = ();
        $matched_crude_flag        = 0;
        $matched_on_ccs_pass_flag  = 0;
        $found_match_flag          = 0;

        $found_unannotated_flag    = 0;
    
        #foreach $ppm_tol (@ppm_tol_array)
        for ($k = @ppm_tol_array - 1; $k < @ppm_tol_array; $k++)
        {
            # be excessively lenient in intial m/z matches
            # we will re-filter them properly later
            $ppm_tol = 2 * $ppm_tol_array[$k];
        
            for ($i = $scils_mz_matched_start; $i < @scils_mz_matched_array; $i++)
            {
                $scils_mz = $scils_mz_matched_array[$i];

                # exit current search, well beyond any sane PPM range
                if ($scils_mz >= $cardinal_mz + 1)
                {
                    last;
                }
                # where to start scanning next time
                if ($scils_mz <= $cardinal_mz - 1)
                {
                    $scils_mz_matched_start = $i;
                }

                $ppm = 1000000 *
                       abs(($cardinal_mz - $scils_mz) / $cardinal_mz);
                
                if ($ppm <= $ppm_tol + 1.0e-5)
                {
                    $crude_mz_match_hash{$scils_mz} = $ppm_tol;
                    $matched_crude_flag = 1;
                }
            }
            
            # found a match, don't try any higher ppm tolerances
            if ($matched_crude_flag)
            {
                last;
            }
        }
        
        # matched to vendor m/z, refine matches with cardinal m/z
        if ($matched_crude_flag)
        {
            @match_scils_mz_array = sort {$a <=> $b} keys
                %crude_mz_match_hash;
            
            # 1st pass, only include CCS matches
            foreach $ppm_tol (@ppm_tol_array)
            {
                foreach $scils_mz (@match_scils_mz_array)
                {
                    $adduct = $scils_ccs_hash{$scils_mz}{adduct};
                    $ccs    = $scils_ccs_hash{$scils_mz}{ccs};

                    # guess missing adduct
                    if (!($adduct =~ /[A-Za-z]/) && $num_assumed_adduct_types)
                    {
                        $adduct = $assume_adduct_array[0];
                    }
    
                    # skip unsupported adducts
                    if (!defined($adduct_offset_hash{$adduct})) { next; }
                    
                    # skip missing CCS
                    if (!defined($ccs) || !is_number($ccs)) { next; }

                    @merged_id_array_1 =
                        keys %{$scils_match_hash{$scils_mz}};

                    @merged_id_array_2 =
                        keys %{$original_lmaps_matches{$scils_mz}};
                    
                    # doesn't matter at this point if we have dupes or not
                    @merged_id_array = sort cmp_merged_id
                        (@merged_id_array_1, @merged_id_array_2);

                    foreach $merged_id (@merged_id_array)
                    {
                        # only include CCS matches this pass
                        $ccs_flag = 0;
                        if (defined($scils_match_hash{$scils_mz}{$merged_id}))
                        {
                            $ccs_flag =
                                $scils_match_hash{$scils_mz}{$merged_id}{matched_on_ccs};
                        }
                        elsif (defined($original_lmaps_matches{$cardinal_mz}) &&
                               defined($original_lmaps_matches{$cardinal_mz}{$merged_id}))
                        {
                            $ccs_flag =
                                $original_lmaps_matches{$cardinal_mz}{$merged_id};
                        }
                        if ($ccs_flag == 0)
                        {
                            next;
                        }


                        $mass_db       = '';
                        $formal_charge = 0;

                        if ($merged_id =~ /^HMDB/)
                        {
                            $mass_db       = $hmdb_ccs_hash{$merged_id}{mass};
                            $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

                            if (!defined($formal_charge) || !is_number($formal_charge))
                            {
                                $formal_charge = 0;
                            }
                        }
                        elsif ($merged_id =~ /^LM/)
                        {
                            $mass_db = $lmaps_ccs_hash{$merged_id}{mass};
                        }


                        if (!is_number($mass_db)) { next; }

                        # skip formal charges, since they are more likely in LC than MALDI
                        if ($formal_charge) { next; }


                        $adduct_offset = $adduct_offset_hash{$adduct};
                        if (!defined($adduct_offset)) { next; }

                        # deal with formal charges
                        if ($formal_charge)
                        {
                            $adduct_charge = $adduct_charge_hash{$adduct};

                            # can never switch between opposite charges, skip it
                            if (($formal_charge > 0 && $adduct_charge < 0) ||
                                ($formal_charge < 0 && $adduct_charge > 0))
                            {
                                next;
                            }
                        
                            # zero out H adduct offsets
                            if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
                            {
                                $adduct_offset = 0;
                            }
                        
                            $delta_charge   = $adduct_charge - $formal_charge;
                            $adduct_offset += $delta_charge * $mass_proton;
                        }


                        $mz_db = $mass_db + $adduct_offset;
                        $ppm   = 1000000 * abs(($mz_db - $cardinal_mz) /
                                               $mz_db);

                        if ($ppm <= $ppm_tol + 1.0e-5)
                        {
                            $cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}{$merged_id} = 1;

                            if ($adduct =~ /[A-Za-z]/)
                            {
                                $cardinal_refined_mz_match_hash{$cardinal_mz}{adduct}{$adduct}   = 1;
                            }
                            if (is_number($ccs))
                            {
                                $cardinal_refined_mz_match_hash{$cardinal_mz}{ccs}{$ccs}         = 1;
                            }

                            $cardinal_refined_mz_match_hash{$cardinal_mz}{ppm_tol}               = $ppm_tol;
                            $cardinal_refined_mz_match_hash{$cardinal_mz}{matched_on_ccs}        = 1;

                            $found_match_flag = 1;
                        }
                    }
                }
                
                if ($found_match_flag)
                {
                    $matched_on_ccs_pass_flag = 1;

                    last;
                }
            }

            # 2nd pass, do not check for CCS matches
            if ($found_match_flag == 0)
            {
              foreach $ppm_tol (@ppm_tol_array)
              {
                foreach $scils_mz (@match_scils_mz_array)
                {
                    $adduct = $scils_ccs_hash{$scils_mz}{adduct};
                    $ccs    = $scils_ccs_hash{$scils_mz}{ccs};

                    # guess missing adduct
                    if (!($adduct =~ /[A-Za-z]/) && $num_assumed_adduct_types)
                    {
                        $adduct = $assume_adduct_array[0];
                    }
    
                    # skip unsupported adducts
                    if (!defined($adduct_offset_hash{$adduct})) { next; }

                    @merged_id_array_1 =
                        keys %{$scils_match_hash{$scils_mz}};

                    @merged_id_array_2 =
                        keys %{$original_lmaps_matches{$scils_mz}};
                    
                    # doesn't matter at this point if we have dupes or not
                    @merged_id_array = sort cmp_merged_id
                        (@merged_id_array_1, @merged_id_array_2);

                    foreach $merged_id (@merged_id_array)
                    {
                        $mass_db       = '';
                        $formal_charge = 0;
                        
                        if ($merged_id =~ /^HMDB/)
                        {
                            $mass_db       = $hmdb_ccs_hash{$merged_id}{mass};
                            $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

                            if (!defined($formal_charge) || !is_number($formal_charge))
                            {
                                $formal_charge = 0;
                            }
                        }
                        elsif ($merged_id =~ /^LM/)
                        {
                            $mass_db = $lmaps_ccs_hash{$merged_id}{mass};
                        }


                        if (!is_number($mass_db)) { next; }

                        # skip formal charges, since they are more likely in LC than MALDI
                        if ($formal_charge) { next; }


                        $adduct_offset = $adduct_offset_hash{$adduct};
                        if (!defined($adduct_offset)) { next; }

                        # deal with formal charges
                        if ($formal_charge)
                        {
                            $adduct_charge = $adduct_charge_hash{$adduct};

                            # can never switch between opposite charges, skip it
                            if (($formal_charge > 0 && $adduct_charge < 0) ||
                                ($formal_charge < 0 && $adduct_charge > 0))
                            {
                                next;
                            }
                        
                            # zero out H adduct offsets
                            if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
                            {
                                $adduct_offset = 0;
                            }
                        
                            $delta_charge   = $adduct_charge - $formal_charge;
                            $adduct_offset += $delta_charge * $mass_proton;
                        }


                        $mz_db = $mass_db + $adduct_offset;
                        $ppm   = 1000000 * abs(($mz_db - $cardinal_mz) /
                                               $mz_db);

                        if ($ppm <= $ppm_tol + 1.0e-5)
                        {
                            $cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}{$merged_id} = 1;

                            if ($adduct =~ /[A-Za-z]/)
                            {
                                $cardinal_refined_mz_match_hash{$cardinal_mz}{adduct}{$adduct}   = 1;
                            }
                            if (is_number($ccs))
                            {
                                $cardinal_refined_mz_match_hash{$cardinal_mz}{ccs}{$ccs}         = 1;
                            }

                            $cardinal_refined_mz_match_hash{$cardinal_mz}{ppm_tol}               = $ppm_tol;
                            $cardinal_refined_mz_match_hash{$cardinal_mz}{matched_on_ccs}        = 0;

                            $found_match_flag = 1;
                        }
                    }
                }
                
                if ($found_match_flag)
                {
                    last;
                }
              }
            }
        }

        if ($found_match_flag)
        {        
            @kept_match_array =
                sort keys %{$cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}};
        }
        # no annotation matches, try matches to unannotated vendor m/z
        else
        {
            foreach $ppm_tol (@ppm_tol_array)
            {
                for ($i = $scils_mz_all_start; $i < @scils_mz_all_array;
                     $i++)
                {
                    $scils_mz = $scils_mz_all_array[$i];
                    $adduct    = $scils_ccs_hash{$scils_mz}{adduct};
                    $ccs       = $scils_ccs_hash{$scils_mz}{ccs};

                    # guess missing adduct
                    if (!($adduct =~ /[A-Za-z]/) && $num_assumed_adduct_types)
                    {
                        $adduct = $assume_adduct_array[0];
                    }
    
                    # skip unsupported adducts
                    if (!defined($adduct_offset_hash{$adduct})) { next; }
                    
                    # skip missing CCS
                    # we assume CCS caused it to not be annotated earlier,
                    # so missing CCS violates that assumption
                    if (!defined($ccs) || !is_number($ccs)) { next; }

                    # exit current search, well beyond any sane PPM range
                    if ($scils_mz >= $cardinal_mz + 1)
                    {
                        last;
                    }
                    # where to start scanning next time
                    if ($scils_mz <= $cardinal_mz - 1)
                    {
                        $scils_mz_all_start = $i;
                    }

                    $ppm = 1000000 *
                           abs(($cardinal_mz - $scils_mz) / $cardinal_mz);

                    if ($ppm <= $ppm_tol + 1.0e-5)
                    {
                        $unannotated_mz_match_hash{$cardinal_mz}{ppm_tol} =
                            $ppm_tol;
                        $unannotated_mz_match_hash{$cardinal_mz}{adduct}{$adduct} =
                            1;

                        # only store CCS if it isn't blank
                        if (is_number($ccs))
                        {
                          $unannotated_mz_match_hash{$cardinal_mz}{ccs}{$ccs} =
                            1;
                        }

                        $found_unannotated_flag = 1;
                    }
                }
                
                if ($found_unannotated_flag)
                {
                    last;
                }
            }
        }

        
        # keep only HMDB identifiers when LipidMaps is the same molecule
        # dedupe on inchi, hmdb, kegg, cid
        if (@kept_match_array)
        {
            %temp_lmaps_dupe_hash  = ();    # is LipidMaps covered by HMDB hit?
            %temp_hmdb_hash        = ();
            %temp_hmdb_inchi_hash  = ();
            %temp_hmdb_kegg_hash   = ();
            %temp_hmdb_cid_hash    = ();
            %temp_lmaps_hash       = ();
            %temp_lmaps_inchi_hash = ();
            %temp_lmaps_kegg_hash  = ();
            %temp_lmaps_cid_hash   = ();
            
            foreach $merged_id (@kept_match_array)
            {
                #printf "%s\t%s\n", $cardinal_mz, $merged_id;

                if ($merged_id =~ /^HMDB/)
                {
                    $temp_hmdb_hash{$merged_id} = 1;
                    
                    $inchi = $hmdb_ccs_hash{$merged_id}{inchi};
                    $kegg  = $hmdb_ccs_hash{$merged_id}{kegg};
                    $cid   = $hmdb_ccs_hash{$merged_id}{cid};
                    
                    if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_inchi_hash{$inchi} = 1;
                    }
                    if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_kegg_hash{$kegg} = 1;
                    }
                    if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_cid_hash{$cid} = 1;
                    }
                }
                elsif ($merged_id =~ /^LM/)
                {
                    $temp_lmaps_hash{$merged_id} = 1;

                    $hmdb  = $lmaps_ccs_hash{$merged_id}{hmdb};
                    $inchi = $lmaps_ccs_hash{$merged_id}{inchi};
                    $kegg  = $lmaps_ccs_hash{$merged_id}{kegg};
                    $cid   = $lmaps_ccs_hash{$merged_id}{cid};
                    
                    if (defined($hmdb) && $hmdb =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_hmdb_hash{$merged_id}{$hmdb} = 1;
                    }
                    if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_inchi_hash{$merged_id}{$inchi} = 1;
                    }
                    if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_kegg_hash{$merged_id}{$kegg} = 1;
                    }
                    if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_cid_hash{$merged_id}{$cid} = 1;
                    }
                }
            }
            
            @temp_hmdb_array       = sort keys %temp_hmdb_hash;
            @temp_lmaps_array      = sort keys %temp_lmaps_hash;
            @temp_lmaps_dupe_array = ();

            foreach $lmaps_id (@temp_lmaps_array)
            {
                $cardinal_lmaps_matches{$cardinal_mz}{$lmaps_id} = 1;
            }

            if (@temp_hmdb_array && @temp_lmaps_array)
            {
                foreach $lmaps_id (@temp_lmaps_array)
                {
                    @temp_array = sort keys %{$temp_lmaps_hmdb_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_inchi_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_inchi_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_kegg_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_kegg_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_cid_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_cid_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }
                }
            }
            
            @temp_lmaps_dupe_array = sort keys %temp_lmaps_dupe_hash;
            foreach $lmaps_id (@temp_lmaps_dupe_array)
            {
                delete
                    $cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}{$lmaps_id};
            }
        }
    }
}


# scan cardinal peaks not present in vendor peaks,
# attempt to annotate them with assumed adduct
@cardinal_mz_novel_array = ();
$total_novel             = 0;
if (defined($cardinal_filename))
{
    @cardinal_mz_array = sort {$a <=> $b} keys %cardinal_seen_mz_hash;
    $num_novel_peaks   = 0;
    
    # store novel cardinal peaks to scan
    foreach $cardinal_mz (@cardinal_mz_array)
    {
        if (!defined($cardinal_refined_mz_match_hash{$cardinal_mz}) &&
            !defined($unannotated_mz_match_hash{$cardinal_mz}))
        {
            $cardinal_mz_novel_array[$num_novel_peaks++] = $cardinal_mz;
        }
    }

    # attempt to annotate novel peaks
    foreach $cardinal_mz (@cardinal_mz_novel_array)
    {
        $mz_floor   = floor($cardinal_mz);
    
        $found_match_flag = 0;
        foreach $ppm_tol (@ppm_tol_array)
        {
            foreach $adduct (@assume_adduct_array)
            {
                # will include MANY outside the acceptable tolerance,
                # but is quick and dirty to code up
                #
                @crude_candidate_array = ();
                if (defined($mz_id_bins_hash{$mz_floor}))
                {
                    @crude_candidate_array =
                        sort keys %{$mz_id_bins_hash{$mz_floor}};
                }
                
                #$count = @crude_candidate_array;
                #printf STDERR "NOVEL_CRUDE\t%s\t%s\t%s\n",
                #    $cardinal_mz, $adduct, $count;

                # refine crude candidate list by m/z tolerance
                foreach $merged_id (@crude_candidate_array)
                {
                    $mass_db       = '';
                    $formal_charge = 0;
                    
                    if ($merged_id =~ /^HMDB/)
                    {
                        $mass_db       = $hmdb_ccs_hash{$merged_id}{mass};
                        $formal_charge = $hmdb_ccs_hash{$merged_id}{charge_formal};

                        if (!defined($formal_charge) || !is_number($formal_charge))
                        {
                            $formal_charge = 0;
                        }
                    }
                    elsif ($merged_id =~ /^LM/)
                    {
                        $mass_db = $lmaps_ccs_hash{$merged_id}{mass};
                    }
                    

                    if (!is_number($mass_db)) { next; }

                    # skip formal charges, since they are more likely in LC than MALDI
                    if ($formal_charge) { next; }


                    $has_ccs_flag = 0;
                    undef $ccs_db;
                    if ($merged_id =~ /^HMDB/)
                    {
                        if (defined($hmdb_ccs_hash{$merged_id}{has_ccs}))
                        {
                            $ccs_db = $hmdb_ccs_hash{$merged_id}{$adduct};

                            $has_ccs_flag = 1;
                        }
                    }
                    elsif ($merged_id =~ /^LM/)
                    {
                        if (defined($lmaps_ccs_hash{$merged_id}{has_ccs}))
                        {
                            $ccs_db = $lmaps_ccs_hash{$merged_id}{$adduct};

                            $has_ccs_flag = 1;
                        }
                    }
                    
                    # skip molecules that have CCS but not for this adduct
                    if ($has_ccs_flag && !defined($ccs_db))
                    {
                        next;
                    }
                    
                    
                    $adduct_offset = $adduct_offset_hash{$adduct};
                    if (!defined($adduct_offset)) { next; }

                    # deal with formal charges
                    if ($formal_charge)
                    {
                        $adduct_charge = $adduct_charge_hash{$adduct};

                        # can never switch between opposite charges, skip it
                        if (($formal_charge > 0 && $adduct_charge < 0) ||
                            ($formal_charge < 0 && $adduct_charge > 0))
                        {
                            next;
                        }
                    
                        # zero out H adduct offsets
                        if ($adduct eq '[M+H]+' || $adduct eq '[M-H]-')
                        {
                            $adduct_offset = 0;
                        }
                    
                        $delta_charge   = $adduct_charge - $formal_charge;
                        $adduct_offset += $delta_charge * $mass_proton;
                    }

                    
                    $mz_db = $mass_db + $adduct_offset;
                    $ppm   = 1000000 * abs(($mz_db - $cardinal_mz) / $mz_db);
                    
                    if ($ppm <= $ppm_tol + 1.0e-5)
                    {
                        $cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}{$merged_id} = 1;

                        if ($adduct =~ /[A-Za-z]/)
                        {
                            $cardinal_refined_mz_match_hash{$cardinal_mz}{adduct}{$adduct}   = 1;
                        }

                        $cardinal_refined_mz_match_hash{$cardinal_mz}{ppm_tol}               = $ppm_tol;
                        $cardinal_refined_mz_match_hash{$cardinal_mz}{matched_on_ccs}        = 0;

                        $found_match_flag = 1;
                    }
                }
            }

            if ($found_match_flag)
            {
                $ppm_tol_kept = $ppm_tol;

                last;
            }
        }
        
        @kept_match_array = ();
        if (defined($cardinal_refined_mz_match_hash{$cardinal_mz}) &&
            defined($cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}))
        {
            @kept_match_array = sort cmp_merged_id keys
                %{$cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}};
        }

        #$count = @kept_match_array;
        #if ($count)
        #{
        #    printf STDERR "NOVEL_UNPRUNED\t%s\t%s\t%s\n",
        #        $cardinal_mz, $ppm_tol_kept, $count;
        #}


        # keep only HMDB identifiers when LipidMaps is the same molecule
        # dedupe on inchi, hmdb, kegg, cid
        if (@kept_match_array)
        {
            %temp_lmaps_dupe_hash  = ();    # is LipidMaps covered by HMDB hit?
            %temp_hmdb_hash        = ();
            %temp_hmdb_inchi_hash  = ();
            %temp_hmdb_kegg_hash   = ();
            %temp_hmdb_cid_hash    = ();
            %temp_lmaps_hash       = ();
            %temp_lmaps_inchi_hash = ();
            %temp_lmaps_kegg_hash  = ();
            %temp_lmaps_cid_hash   = ();

            foreach $merged_id (@kept_match_array)
            {
                if ($merged_id =~ /^HMDB/)
                {
                    $temp_hmdb_hash{$merged_id} = 1;
                    
                    $inchi = $hmdb_ccs_hash{$merged_id}{inchi};
                    $kegg  = $hmdb_ccs_hash{$merged_id}{kegg};
                    $cid   = $hmdb_ccs_hash{$merged_id}{cid};
                    
                    if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_inchi_hash{$inchi} = 1;
                    }
                    if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_kegg_hash{$kegg} = 1;
                    }
                    if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                    {
                        $temp_hmdb_cid_hash{$cid} = 1;
                    }
                }
                elsif ($merged_id =~ /^LM/)
                {
                    $temp_lmaps_hash{$merged_id} = 1;

                    $hmdb  = $lmaps_ccs_hash{$merged_id}{hmdb};
                    $inchi = $lmaps_ccs_hash{$merged_id}{inchi};
                    $kegg  = $lmaps_ccs_hash{$merged_id}{kegg};
                    $cid   = $lmaps_ccs_hash{$merged_id}{cid};
                    
                    if (defined($hmdb) && $hmdb =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_hmdb_hash{$merged_id}{$hmdb} = 1;
                    }
                    if (defined($inchi) && $inchi =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_inchi_hash{$merged_id}{$inchi} = 1;
                    }
                    if (defined($kegg) && $kegg =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_kegg_hash{$merged_id}{$kegg} = 1;
                    }
                    if (defined($cid) && $cid =~ /[A-Za-z0-9]/)
                    {
                        $temp_lmaps_cid_hash{$merged_id}{$cid} = 1;
                    }
                }
            }
            
            @temp_hmdb_array       = sort keys %temp_hmdb_hash;
            @temp_lmaps_array      = sort keys %temp_lmaps_hash;
            @temp_lmaps_dupe_array = ();

            foreach $lmaps_id (@temp_lmaps_array)
            {
                $cardinal_lmaps_matches{$cardinal_mz}{$lmaps_id} = 0;
            }

            if (@temp_hmdb_array && @temp_lmaps_array)
            {
                foreach $lmaps_id (@temp_lmaps_array)
                {
                    @temp_array = sort keys %{$temp_lmaps_hmdb_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_inchi_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_inchi_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_kegg_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_kegg_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }

                    @temp_array = sort keys %{$temp_lmaps_cid_hash{$lmaps_id}};
                    foreach $identifier (@temp_array)
                    {
                        if (defined($temp_hmdb_cid_hash{$identifier}))
                        {
                            $temp_lmaps_dupe_hash{$lmaps_id} = 1;
                            next;
                        }
                    }
                }
            }
            
            @temp_lmaps_dupe_array = sort keys %temp_lmaps_dupe_hash;
            foreach $lmaps_id (@temp_lmaps_dupe_array)
            {
                delete
                    $cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}{$lmaps_id};
            }
        }
        
        # re-build kept match array after deleting LipidMaps HMDB dupes
        @kept_match_array = sort
            keys %{$cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}};
        $count            = @kept_match_array;
        
        if ($count > 0)
        {
            $total_novel++;

            #printf STDERR "NOVEL\t%.4f\t%s\t%s\n",
            #    $cardinal_mz, $ppm_tol_kept, $count;
        }
    }
}

# output cardinal m/z table
if (defined($cardinal_filename))
{
  if ($total_novel)
  {
      printf STDERR "Annotated novel Cardinal peak count:\t%d\n", $total_novel;
  }


  @array = split /\t/, $cardinal_header_line;
  
  for ($i_row = 0; $i_row < @array; $i_row++)
  {
      $field = $array[$i_row];
  
      # insert new annotation columns
      if ($i_row == $cardinal_mz_col)
      {
          if ($i_row)
          {
              printf "\t";
          }

          printf "%s",   'm/z';
          printf "\t%s", 'row order original';
          printf "\t%s", 'm/z identifier';
          printf "\t%s", 'Adduct';
          printf "\t%s", 'CCS';
          printf "\t%s", 'PPM tolerence';
          printf "\t%s", 'CCS matched';
          printf "\t%s", 'Number of matches';
          printf "\t%s", 'Formula';
          printf "\t%s", 'Name';
          printf "\t%s", 'InChI';
          printf "\t%s", 'HMDB';
          printf "\t%s", 'KEGG';
          printf "\t%s", 'PubChemCID';
          printf "\t%s", 'LipidMaps';    # will include removed HMDB dupes

          # insert additional annotation
          for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
          {
              $col = $more_hmdb_annotation_cols[$i];
              
              if (defined($col))
              {
                  printf "\t%s",
                      unicode_to_ascii($more_hmdb_annotation_headers[$i]);
              }
          }
          for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
          {
              $col = $more_lmaps_annotation_cols[$i];
              
              if (defined($col))
              {
                  printf "\t%s",
                      unicode_to_ascii($more_lmaps_annotation_headers[$i]);
              }
          }
      }
      else
      {
          printf "\t%s", $field;
      }
  }
  printf "\n";


  $row_output = 0;
  for ($row = 0; $row < @cardinal_row_array; $row++)
  {
    %seen_formula_hash        = ();
    %seen_name_conformed_hash = ();
    %seen_inchi_hash          = ();
    %seen_hmdb_hash           = ();
    %seen_kegg_hash           = ();
    %seen_cid_hash            = ();
    
    %seen_more_hmdb_hash      = ();
    %seen_more_lmaps_hash     = ();

    $ccs_str     = '';
    $adduct_str  = '';
    $formula_str = '';
    $name_str    = '';
    $inchi_str   = '';
    $hmdb_str    = '';
    $kegg_str    = '';
    $cid_str     = '';
    $lmaps_str   = '';
    
    # initialize more annotation strings
    for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
    {
        $more_hmdb_str_array[$i]  = '';
    }
    for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
    {
        $more_lmaps_str_array[$i] = '';
    }

    
    $num_matches    = 0;
    $ppm_tol        = 0;
    $matched_on_ccs = 0;


    $cardinal_mz = $cardinal_row_array[$row]{mz};
    $line        = $cardinal_row_array[$row]{line};
    @array       = split /\t/, $line, -1;


    if (is_number($cardinal_mz) &&
        defined($cardinal_refined_mz_match_hash{$cardinal_mz}))
    {
        if ($cardinal_refined_mz_match_hash{$cardinal_mz}{ppm_tol})
        {
            $ppm_tol = $cardinal_refined_mz_match_hash{$cardinal_mz}{ppm_tol};
        }
        if ($cardinal_refined_mz_match_hash{$cardinal_mz}{matched_on_ccs})
        {
            $matched_on_ccs = 1;
        }
        
        @temp_array = sort 
            keys %{$cardinal_refined_mz_match_hash{$cardinal_mz}{adduct}};
        $adduct_str = join ' | ', @temp_array;

        @temp_array = sort {$a <=> $b}
            keys %{$cardinal_refined_mz_match_hash{$cardinal_mz}{ccs}};
        $ccs_str = join ' | ', @temp_array;


        @merged_id_array = sort cmp_merged_id
            keys %{$cardinal_refined_mz_match_hash{$cardinal_mz}{merged_id}};
        
        $num_matches = @merged_id_array;
        
        for ($i = 0; $i < $num_matches && $i < $max_hits; $i++)
        {
            $merged_id = $merged_id_array[$i];

            $formula = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $formula = $hmdb_ccs_hash{$merged_id}{formula};
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $formula = $lmaps_ccs_hash{$merged_id}{formula};
            }
            if ($formula =~ /[A-Za-z0-9]/ &&
                !defined($seen_formula_hash{$formula}))
            {
                if ($formula_str eq '')
                {
                    $formula_str = $formula;
                }
                else
                {
                    $formula_str .= ' | ' . $formula;
                }
                $seen_formula_hash{$formula} = 1;
            }

            $name = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $name = $hmdb_ccs_hash{$merged_id}{name};
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $name = $lmaps_ccs_hash{$merged_id}{name};
            }
            $name_conformed = conform_name($name);
            if ($name =~ /[A-Za-z0-9]/ &&
                !defined($seen_name_conformed_hash{$name_conformed}))
            {
                if ($name_str eq '')
                {
                    $name_str = $name;
                }
                else
                {
                    $name_str .= ' | ' . $name;
                }
                $seen_name_conformed_hash{$name_conformed} = 1;
            }

            $inchi = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $inchi = $hmdb_ccs_hash{$merged_id}{inchi};
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $inchi = $lmaps_ccs_hash{$merged_id}{inchi};
            }
            if ($inchi =~ /[A-Za-z0-9]/ &&
                !defined($seen_inchi_hash{$inchi}))
            {
                if ($inchi_str eq '')
                {
                    $inchi_str = $inchi;
                }
                else
                {
                    $inchi_str .= ' | ' . $inchi;
                }
                $seen_inchi_hash{$inchi} = 1;
            }

            $hmdb = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $hmdb = $merged_id;
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $hmdb = $lmaps_ccs_hash{$merged_id}{hmdb};
            }
            if ($hmdb =~ /[A-Za-z0-9]/ &&
                !defined($seen_hmdb_hash{$hmdb}))
            {
                if ($hmdb_str eq '')
                {
                    $hmdb_str = $hmdb;
                }
                else
                {
                    $hmdb_str .= ' | ' . $hmdb;
                }
                $seen_hmdb_hash{$hmdb} = 1;
            }

            $kegg = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $kegg = $hmdb_ccs_hash{$merged_id}{kegg};
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $kegg = $lmaps_ccs_hash{$merged_id}{kegg};
            }
            if ($kegg =~ /[A-Za-z0-9]/ &&
                !defined($seen_kegg_hash{$kegg}))
            {
                if ($kegg_str eq '')
                {
                    $kegg_str = $kegg;
                }
                else
                {
                    $kegg_str .= ' | ' . $kegg;
                }
                $seen_kegg_hash{$kegg} = 1;
            }

            $cid = '';
            if ($merged_id =~ /^HMDB/i)
            {
                $cid = $hmdb_ccs_hash{$merged_id}{cid};
            }
            elsif ($merged_id =~ /^LM/i)
            {
                $cid = $lmaps_ccs_hash{$merged_id}{cid};
            }
            if ($cid =~ /[A-Za-z0-9]/ &&
                !defined($seen_cid_hash{$cid}))
            {
                if ($cid_str eq '')
                {
                    $cid_str = $cid;
                }
                else
                {
                    $cid_str .= ' | ' . $cid;
                }
                $seen_cid_hash{$cid} = 1;
            }
        }

        @lmaps_id_array = ();
        if (defined($cardinal_lmaps_matches{$cardinal_mz}))
        {
            @lmaps_id_array = sort cmp_merged_id keys %{$cardinal_lmaps_matches{$cardinal_mz}};
            $lmaps_str      = join ' | ', @lmaps_id_array;
        }

        # more annotations
        foreach $merged_id (@merged_id_array)
        {
            if (!defined($more_hmdb_annotation_hash{$merged_id}))
            {
                next;
            }

            for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
            {
                $value = $more_hmdb_annotation_hash{$merged_id}[$i];

                if (defined($value) && $value =~ /[A-Za-z0-9]/ &&
                    !(defined($seen_more_hmdb_hash{$i}) &&
                      defined($seen_more_hmdb_hash{$i}{$value})))
                {
                    if ($more_hmdb_str_array[$i] eq '')
                    {
                        $more_hmdb_str_array[$i] = $value;
                    }
                    else
                    {
                        $more_hmdb_str_array[$i] .= ' | ' . $value;
                    }

                    $seen_more_hmdb_hash{$i}{$value} = 1;
                }
            }
        }

        foreach $merged_id (@lmaps_id_array)
        {
            if (!defined($more_lmaps_annotation_hash{$merged_id}))
            {
                next;
            }

            for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
            {
                $value = $more_lmaps_annotation_hash{$merged_id}[$i];

                if (defined($value) && $value =~ /[A-Za-z0-9]/ &&
                    !(defined($seen_more_lmaps_hash{$i}) &&
                      defined($seen_more_lmaps_hash{$i}{$value})))
                {
                    if ($more_lmaps_str_array[$i] eq '')
                    {
                        $more_lmaps_str_array[$i] = $value;
                    }
                    else
                    {
                        $more_lmaps_str_array[$i] .= ' | ' . $value;
                    }

                    $seen_more_lmaps_hash{$i}{$value} = 1;
                }
            }
        }
    }
    
    # blank out zero matches
    if ($num_matches == 0)
    {
        $num_matches    = '';
        $ccs_str        = '';
        $ppm_tol        = '';
        $matched_on_ccs = '';
    }

    # no annotation hits, check for unannotated vendor mz hits
    if (is_number($cardinal_mz) &&
        defined($unannotated_mz_match_hash{$cardinal_mz}))
    {
        if ($unannotated_mz_match_hash{$cardinal_mz}{ppm_tol})
        {
            $ppm_tol = $unannotated_mz_match_hash{$cardinal_mz}{ppm_tol};
        }
        
        @temp_array = sort 
            keys %{$unannotated_mz_match_hash{$cardinal_mz}{adduct}};
        $adduct_str = join ' | ', @temp_array;

        # CCS can be blank, so sort alphanumerically
        $ccs_str = '';
        @temp_array = sort { cmp_args_alphanumeric($a, $b) }
            keys %{$unannotated_mz_match_hash{$cardinal_mz}{ccs}};
        if (@temp_array)
        {
            $ccs_str = join ' | ', @temp_array;

            # strip blank CCS from ends
            $ccs_str =~ s/^\s*\|\s*+//;
            $ccs_str =~ s/\s*\|+\s*$//;
        }
    }


    $row_output++;

    for ($i_row = 0; $i_row < @array; $i_row++)
    {
        $field = $array[$i_row];
  
        # insert new annotation columns
        if ($i_row == $cardinal_mz_col)
        {
            if ($i_row)
            {
                printf "\t";
            }

            printf "%s",   $cardinal_mz;
            printf "\t%s", $row_output;
            printf "\t%s", $array[$cardinal_mz_col];
            printf "\t%s", $adduct_str;
            printf "\t%s", $ccs_str;
            printf "\t%s", $ppm_tol;
            printf "\t%s", $matched_on_ccs;
            printf "\t%s", $num_matches;
            printf "\t%s", $formula_str;
            printf "\t%s", unicode_to_ascii($name_str);
            printf "\t%s", $inchi_str;
            printf "\t%s", $hmdb_str;
            printf "\t%s", $kegg_str;
            printf "\t%s", $cid_str;
            printf "\t%s", $lmaps_str;

            # insert additional annotation
            for ($i = 0; $i < @more_hmdb_annotation_headers; $i++)
            {
                $col = $more_hmdb_annotation_cols[$i];
                
                if (defined($col))
                {
                    printf "\t%s", unicode_to_ascii($more_hmdb_str_array[$i]);
                }
            }
            for ($i = 0; $i < @more_lmaps_annotation_headers; $i++)
            {
                $col = $more_lmaps_annotation_cols[$i];
                
                if (defined($col))
                {
                    printf "\t%s", unicode_to_ascii($more_lmaps_str_array[$i]);
                }
            }
        }
        else
        {
            printf "\t%s", $field;
        }
    }
    printf "\n";
  }
}
