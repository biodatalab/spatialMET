#!/usr/bin/perl -w


# HMDB names with corrupt unicode characters:
# HMDB0300900    corrupt name: (2E_4Z)????decadienoyl-CoA
#                http://qpmf.rx.umaryland.edu/PAMDB?MetID=PAMDB001410
#                should probably be (2E_4Z)-decadienoyl-CoA
# HMDB0251069    corrupt name: 2,2???-(Hydroxynitrosohydrazino)bis-ethanamine
#                bloodexposome.org: 2,2'-( ...
# HMDB0250632    corrupt name: ... cyclic (3?5)-disulfide
#                bloodexposome.org: ... cyclic (35)-disulfide
# HMDB0304570    corrupt name: 4?-O-glucoside;   substitute with '
# HMDB0304569    corrupt name: 2??-O-rhamnoside; substitute with '
# HMDB0304547    corrupt name: ?-hydroxylaurate; substitute with omega
# HMDB0242122    corrupt traditional iupac; substitute <?> char with alpha

# HMDB names that are just plain nonsense
# we'll want to use the traditional IUPAC instead
# HMDB0302501    "5"


use Encode qw(decode encode);
# !!!!!!!!!!!!!!!!!!!
# !!! Dependency: !!!
# !!!             !!!
# !!!   https://metacpan.org/dist/Text-Unidecode
# !!!             !!!
# !!!!!!!!!!!!!!!!!!!
use Text::Unidecode;

sub cmp_unicode_hex_str
{
    my $len_a = length $a;
    my $len_b = length $b;
    
    if ($len_a < $len_b) { return -1; }
    if ($len_a > $len_b) { return  1; }
    
    return $a cmp $b;
}


$header_keep_order_hash{'accession'}                       = 0;
$header_keep_order_hash{'predicted_properties::mono_mass'} = 1;
$header_keep_order_hash{'chemical_formula'}                = 2;
$header_keep_order_hash{'predicted_properties::formal_charge'} = 3;
$header_keep_order_hash{'predicted_properties::physiological_charge'} = 4;
$header_keep_order_hash{'smiles'}                          = 5;
$header_keep_order_hash{'name'}                            = 10;
$header_keep_order_hash{'traditional_iupac'}               = 11;
$header_keep_order_hash{'iupac_name'}                      = 12;
$header_keep_order_hash{'synonyms::synonym'}               = 13;
$header_keep_order_hash{'inchi'}                           = 50;
$header_keep_order_hash{'inchikey'}                        = 51;
#$header_keep_order_hash{'bigg_id'}                         = 100;
#$header_keep_order_hash{'biocyc_id'}                       = 101;
$header_keep_order_hash{'cas_registry_number'}             = 102;
$header_keep_order_hash{'chebi_id'}                        = 103;
#$header_keep_order_hash{'chemspider_id'}                   = 104;
$header_keep_order_hash{'drugbank_id'}                     = 105;
#$header_keep_order_hash{'foodb_id'}                        = 106;
$header_keep_order_hash{'kegg_id'}                         = 107;
$header_keep_order_hash{'biological_properties::pathways::pathway::kegg_map_id'} = 107.1;
#$header_keep_order_hash{'knapsack_id'}                     = 108;
#$header_keep_order_hash{'metlin_id'}                       = 109;
$header_keep_order_hash{'pubchem_compound_id'}             = 110;
#$header_keep_order_hash{'vmh_id'}                          = 111;

$header_keep_order_hash{'taxonomy::super_class'}           = 200;
$header_keep_order_hash{'taxonomy::class'}                 = 201;
$header_keep_order_hash{'taxonomy::sub_class'}             = 202;
$header_keep_order_hash{'taxonomy::direct_parent'}         = 203;
$header_keep_order_hash{'taxonomy::molecular_framework'}   = 204;


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
#$unicode_to_ascii_hash{"\x{FFFD}"} = 'a';


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

            # bogus character, probably supposed to be alpha?
            # it only occurs once in all of the HMDB compound name fields
            # HMDB0242122 ?-D-galactopyranoside, ethyl
            #
            # restrict alpha substitution to only HMDB0242122,
            # just in case more appear elsewhere in the future
            #
            if (!defined($c_new) &&
                $c_single  eq "\x{FFFD}" &&
                $accession eq 'HMDB0242122')
            {
                $c_new = 'a';
            }

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


# lookup table for converting numbers into ordered greek letters
$number_letter_hash{1} = 'a';
$number_letter_hash{2} = 'b';
$number_letter_hash{3} = 'g';
$number_letter_hash{4} = 'd';
$number_letter_hash{5} = 'e';
$number_letter_hash{6} = 'z';
$number_letter_hash{7} = 'h';

sub preprocess_name
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

    # replace Greek letters at letter boundaries
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
    $name =~ s/acid,(.*)ic$/$1ic acid/;   # reorder weird MeSH, HMDB notation
    $name =~ s/anoic/yric/g;              # Butanoic acid --> Butyric acid
    $name =~ s/anoate\b/yrate/g;          # Butanoate     --> Butyrate
    $name =~ s/ic acid\b/ate/g;           # Glutamic acid --> Glutamate
    
    # strip the L- from L-aminoacids
    # only on word boundary, so we don't strip DL-aminoacid
    $name =~ s/\bl-(alanine)/$1/g;
    $name =~ s/\bl-(arginine)/$1/g;
    $name =~ s/\bl-(asparagine)/$1/g;
    $name =~ s/\bl-(aspartic acid)/$1/g;
    $name =~ s/\bl-(cysteine)/$1/g;
    $name =~ s/\bl-(glutamine)/$1/g;
    $name =~ s/\bl-(glutamic acid)/$1/g;
    # $name =~ s/\bl-(glycine)/$1/g;           # L- is never used !!
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
    $name =~ s/\bsulfid\b/sulfide/g;
    $name =~ s/\bsulphid\b/sulfide/g;
    $name =~ s/\bsulphide\b/sulfide/g;

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
    $name =~ s/mono(\S)/$1/g;
    
    # condense everything that isn't a letter, number, comma, or plus/minus
    # except when between two numbers
    #
    # protect -) as in (+/-) or (-) using capital letters
    $name =~ s/\(-|-\)/MINUS/g;       # protect minus signs
    $name =~ s/[^A-Za-z0-9,+]/-/g;    # convert to hyphens
    $name =~ s/-+/-/g;                # condense multiple hyphens in a row
    $name =~ s/(^-|-$)//g;            # strip leading/trailing hyphens
    $name =~ s/([,+])-/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/-([,+])/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/([a-z])-/$1/g;         # strip hyphens next to letters
    $name =~ s/-([a-z])/$1/g;         # strip hyphens near to letters
    $name =~ s/,+/,/g;                # condense multiple commas in a row
    $name =~ s/\++/\+/g;              # condense multiple pluses in a row
    
    # de-protect and condense minus signs
    #
    # hypothetical example: 1-(-)-galacturonate   -->   1-galacturonate
    #                       1-(+)-galacturonate   -->   1+galacturonate
    #                       1-galacturonate       -->   1galacturonate
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



$hmdb_parsed_filename = shift;

open HMDB, "$hmdb_parsed_filename" or die "ABORT -- cannot open file $hmdb_parsed_filename\n";

$line = <HMDB>;
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;	# we don't expect any real "" in the headers
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;

    if ($array[$i] =~ /\S/)
    {    
        $hmdb_header_col_hash{$array[$i]} = $i;
        # $header_col_array[$i] = $array[$i];
    }
}

$hmdb_accession_col = $hmdb_header_col_hash{'accession'};
@header_array = sort {$header_keep_order_hash{$a} <=>
                      $header_keep_order_hash{$b}}
                      keys %header_keep_order_hash;

# assign column numbers to headers
@header_keep_hash = ();
foreach $header (@header_array)
{
    $header_keep_hash{$header} = $hmdb_header_col_hash{$header};
}
@header_array = sort {$header_keep_order_hash{$a} <=>
                      $header_keep_order_hash{$b}}
                      keys %header_keep_hash;


# read in the rest of the file
while(defined($line=<HMDB>))
{
    $line  =~ s/[\r\n]+//;
    @array = split /\t/, $line, -1;

    for ($i = 0; $i < @array; $i++)
    {
        # continue stripping problematic stuff until all has been stripped
        do
        {
            $changed_flag = 0;
            
            # remove pre-existing escapes or start/end double quotes,
            # since either messes up ="" escapes
            while ($array[$i] =~ s/^\=*\"(.*?)\"$/$1/)
            {
                $changed_flag = 1;
            }
        
            # remove leading ", since they mess up Excel in general
            #
            # this must be done after "", but before leading/trailing spaces,
            # since removing leading/trailing spaces could result in more
            # full "" enclosures, which would then be messed up by removing
            # only the leading "
            #
            while ($array[$i] =~ s/^\"//)
            {
                $changed_flag = 1;
            }

            # remove leading spaces, since they won't protect long numbers,
            # and will cause various REGEX to fail
            if ($array[$i] =~ s/^\s+//)
            {
                $changed_flag = 1;
            }

            # remove trailing spaces, since they won't protect dates,
            # and will cause various REGEX to fail
            if ($array[$i] =~ s/\s+$//)
            {
                $changed_flag = 1;
            }
        } while ($changed_flag);
    }
    
    $accession = $array[$hmdb_accession_col];

    foreach $header (@header_array)
    {
        $col        = $header_keep_hash{$header};
        $value_str = $array[$col];

        # blank fields at end of line weren't output, so line ends early
        if (!defined($value_str))
        {
            next;
        }

        @split_array = split /\|\|/, $value_str;
        
        foreach $value (@split_array)
        {
            if ($value =~ /[A-Za-z0-9]/)
            {
                # HMDB0303339
                #
                # Umm... no.  Lunatone, the metabolite, is not a Pokemon.
                if ($value =~ /List of generation III Pok/i)
                {
                    next;
                }


                # HACK -- deal with brand junk
                $brand = '';
                $trimmed = $value;
                if ($trimmed =~ s/^(\S+)\s+(\S+.*?)\s*\bbrand$/$1/i)
                {
                    $brand = lc $2;
                    $value = $trimmed;
                }
                if ($trimmed =~ s/(.*?)\s*\bbrand(\s*\d+)*\s*(\bof\b)*\s*//i)
                {
                    $brand = lc $1;
                    $value = $trimmed;
                }


                # replace common unicode with ASCII
                $full_unicode = $value;
                $value = unicode_to_ascii($value);

            
                # scan data string for remaining unicode characters
                $unicode = $value;
                $unicode =~ s/[\x00-\x7F]//g;

                
                # HACK -- if we have any left, they are usually corrupt
                # restore all unicode characters so that we can maybe attempt
                # to uncorrupt it later
                if ($unicode ne '')
                {
                    $value  = $full_unicode;
                    $unicode = $value;
                    $unicode =~ s/[\x00-\x7F]//g;
                }

                if ($unicode ne '')
                {
                    # first, decode the unicode string
                    # into single characters, so substr works correctly
                    utf8::decode($unicode);
                    $len = length $unicode;
                    
                    $seen_unicode_value_hash{$value} = 1;

                    for ($j = 0; $j < $len; $j++)
                    {
                        $c_single = substr $unicode, $j, 1;
                    
                        # encode it back again, since input is multi-byte chars
                        $c_multi = $c_single;
                        utf8::encode($c_multi);
                        
                        $utf16_hex_str = unpack("H*", encode("UTF-16BE", $c_single));

                        # convert to hexidecimal for storage and readability
                        $utf8_hex_str = unpack("H*", $c_multi);
                        
                        $utf8_hex_to_utf8_multi_hash{$utf8_hex_str} = $c_multi;
                        $utf8_hex_to_ascii_hash{$utf8_hex_str}      = unidecode($c_single);
                        
                        $utf8_hex_to_utf16_hex_hash{$utf8_hex_str}  =
                            'U+' . uc $utf16_hex_str;
                            
                        # keep track of where we saw it
                        $utf8_hex_seen_header_hash{$utf8_hex_str}{$header} += 1;
                        $seen_unicode_header_hash{$header} += 1;
                        
                        $seen_unicode_header_hash{$header} = 1;
                    }
                }

                # $value = preprocess_name($value);
            
                $hmdb_data_hash{$accession}{$header}{$value} = 1;
                $seen_hmdb_header_hash{$header} = 1;
            }
        }
    }
    
    # given name is nonsense (just a single number by itself)
    #  example: HMDB0302501
    #
    # replace with traditional IUPAC
    @name_array = sort keys %{$hmdb_data_hash{$accession}{name}};
    $name       = $name_array[0];
    if ($name =~ /^[0-9]+$/)
    {
        @temp_array = sort keys
                      %{$hmdb_data_hash{$accession}{traditional_iupac}};

        if (@temp_array && $temp_array[0] =~ /[A-Za-z]/)
        {
            $name_new = $temp_array[0];

            delete $hmdb_data_hash{$accession}{name}{$name};
            $hmdb_data_hash{$accession}{name}{$name_new} = 1;

            printf STDERR "Replacing bogus name with IUPAC traditional:\t%s\t%s\t%s",
                $accession, $name, $name_new;
        }
        else
        {
            @temp_array = sort keys
                          %{$hmdb_data_hash{$accession}{iupac_name}};

            if (@temp_array && $temp_array[0] =~ /[A-Za-z]/)
            {
                $name_new = $temp_array[0];

                delete $hmdb_data_hash{$accession}{name}{$name};
                $hmdb_data_hash{$accession}{name}{$name_new} = 1;

                printf STDERR "Replacing bogus name with IUPAC:\t%s\t%s\t%s",
                    $accession, $name, $name_new;
            }
        }
    }
}
close HMDB;


# print all remaining values containing unicode characters
open UNICODE, ">detected_text_unicode_remaining.txt" or die "can't open unicode output file\n";
foreach $value (sort keys %seen_unicode_value_hash)
{
    print UNICODE "$value\r\n";
}
close UNICODE;


# print list of encountered unicode characters
@seen_unicode_header_array = sort keys %seen_unicode_header_hash;
open UNICODE, ">detected_header_unicode_characters.csv" or die "can't open unicode output file\n";

# output UTF8 byte order mark (BOM),
# which is evil, but the only way to make Excel parse the file correctly
if (@seen_unicode_header_array)
{
    print UNICODE "\xEF\xBB\xBF";
    print UNICODE "\n";
    print UNICODE "# Don't forget to remove the UTF8 byte order mark from the start of this file!\n";

    # tab-delimited UTF8 import is broken in Excel (it removes the tabs !!)
    # so, we have to comma delimit everything

    printf UNICODE "%s,%s,%s,%s",
        'CodePoint', 'UTF-8', 'Unicode', 'ASCII';
    foreach $header (@seen_unicode_header_array)
    {
        # escape commas and double quotes
        if ($header =~ /\"/ || $header =~ /,/)
        {
            $header =~ s/\"/\"\"/g;
            $header = "\"" . $header . "\"";
        }

        printf UNICODE ",%s", $header;
    }
    printf UNICODE "\n";

    foreach $utf8_hex (sort cmp_unicode_hex_str keys %utf8_hex_to_utf8_multi_hash)
    {
        printf UNICODE "%s",  $utf8_hex_to_utf16_hex_hash{$utf8_hex};
        printf UNICODE ",%s", $utf8_hex;
        printf UNICODE ",%s", $utf8_hex_to_utf8_multi_hash{$utf8_hex};

        $ascii = $utf8_hex_to_ascii_hash{$utf8_hex};

        # escape any ASCII newline characters with backslashes
        $ascii =~ s/\r/\\r/g;
        $ascii =~ s/\n/\\n/g;

        # escape commas and double quotes
        if ($ascii =~ /\"/ || $ascii =~ /,/)
        {
            $ascii =~ s/\"/\"\"/g;
            $ascii = "\"" . $ascii . "\"";
        }
        printf UNICODE ",%s", $ascii;

        # print the counts of where we saw it
        foreach $header (@seen_unicode_header_array)
        {
            $count = $utf8_hex_seen_header_hash{$utf8_hex}{$header};
        
            if (!defined($count) || $count == 0)
            {
                $count = "";
            }

            printf UNICODE ",%s", $count;
        }

        printf UNICODE "\n";
    }
}
close UNICODE;


@hmdb_accession_array = sort keys %hmdb_data_hash;
@seen_hmdb_header_array = sort {$header_keep_order_hash{$a} <=>
                                $header_keep_order_hash{$b}}
                          keys %seen_hmdb_header_hash;

$header_line = join "\t", @seen_hmdb_header_array;
$header_line =~ s/[:\w]+:://g;

# print header line
print "$header_line\n";

foreach $accession (@hmdb_accession_array)
{
    for ($i = 0; $i < @seen_hmdb_header_array; $i++)
    {
        $header = $seen_hmdb_header_array[$i];
    
        @value_array = sort keys %{$hmdb_data_hash{$accession}{$header}};
        
        $value = join '|', @value_array;
        
        if ($i)
        {
            print "\t";
        }
        
        print $value;
    }
    print "\n";
}
