=head1 CONTACT

 author: Mariusz Butkiewicz <mariusz.butkiewicz@case.edu>

=cut

=head1 NAME

=head1 SYNOPSIS 

=head1 DESCRIPTION

=cut

package cocos;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::Perl;
use List::MoreUtils qw( pairwise );
use List::MoreUtils qw( each_array );

my $termination_code = 0;

sub version {
    return '84';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        TRANS_RET_PCT => "retained unaltered transcript in percent (COCOS Plugin)",
        TRANS_LEN => "length of unaltered transcript in nucleotide bases (COCOS Plugin)",
        ALT_AA_SEQ => "length of alterated amino acid sequence (COCOS Plugin)",
        TERM_TYPE => "termination type of sequence alteration determination (COCOS Plugin)"       
    };
}

sub run {
    my ($self, $transcript_variation_allele, $vep_line) = @_;

    my $trans = $transcript_variation_allele->transcript;
    my $variant_cdna_pos = $vep_line->{'cDNA_position'};
    my $variant_cds_pos = $vep_line->{'CDS_position'};
    my $transcript_id = $vep_line->{'Feature'};
    my $uploaded_variation = $vep_line->{'Uploaded_variation'};
    my $consequence = $vep_line->{'Consequence'};
    $uploaded_variation =~ s/\//_/g;

    return {} if $consequence !~ /stop_lost/ && $consequence !~ /frame/ || $consequence =~ /NMD/;

    if($variant_cdna_pos !~ /^-$/ && $variant_cds_pos !~ /^-$/) {
        $variant_cdna_pos = (split /-/, $variant_cdna_pos)[0]; 

        my($org_trans_seq, $result) = process_transcript( $transcript_variation_allele);
        my $result_aa = translate_seq_string($result);

	# get chr and pos for current variant
        my $chr = $transcript_variation_allele->transcript->feature_Slice->seq_region_name;
        my $pos = ($transcript_variation_allele->transcript->get_TranscriptMapper->cdna2genomic($variant_cdna_pos,$variant_cdna_pos))[0]->start;
	my $chr_pos = $chr."_".$pos;
        my $allele_str = $transcript_variation_allele->allele_string;
	my $delim = "|";

	my $fasta_header = "chr=".$chr.$delim;
	$fasta_header .= "pos=".$pos.$delim;
        $fasta_header .= "".$transcript_id.$delim;
	$fasta_header .= "".$consequence.$delim;
	$fasta_header .= "var_pos=".(length $org_trans_seq).$delim;
	$fasta_header .= "tr_len=".((length transcript_coding_seq($trans))-3).$delim;
	$fasta_header .= "num_aa=".(length $result_aa).$delim;
	$fasta_header .= "".$allele_str;
	
        write_to_fasta('./', $fasta_header, $result);
        write_to_fasta('./', $fasta_header, $result_aa);

	my $trans = $transcript_variation_allele->transcript;
        my $trans_len_pct_retained = ((length $org_trans_seq) + (length $result)) / (length transcript_coding_seq($trans));	
	    
        return {
            TRANS_LEN => length transcript_coding_seq($trans),
            TRANS_RET_PCT => $trans_len_pct_retained,
            ALT_AA_SEQ => $result_aa,
	    TERM_TYPE => $termination_code,
	    FILE_NAME => $chr_pos
        };
    }

    return {};
}

#write header and seq to fasta file
sub write_to_fasta {
    my ($file_path, $fasta_header, $seq) = @_;
    my $seq_len = length $seq;

    # default fasta length of 80 characters per line
    my $len = 80;
    my $content = ">$fasta_header\n";
    while (my $chunk = substr($seq, 0, $len, "")) {
    	$content .= "$chunk\n";
    }
 
    my $filename = $file_path . "cocos.fasta";
    my $die_msg = "Could not open file '$filename'!";
    $die_msg .= " Please check if file premissions are correct or the appropriate directory exists!";
    open(my $fh, '>>:encoding(UTF-8)', $filename) or die $die_msg;
    print $fh $content;
    close $fh;
}

# translate a DNA sequence string into a amino acid sequence string
# if sequence is not empty or contains 'X' characters (unknown)
sub translate_seq_string {
    my ($custom_seq) = @_;

    if($custom_seq =~ /X/) {
	print "sequence can not be translated! contains X characters!\n";
	print "$custom_seq\n";
        return "";
    }

    if($custom_seq ne '') {
        my $seq_object = new_sequence($custom_seq,'','');
        return $seq_object->translate->seq();
    }

    return '';
}

#get coding seq excluding 5'UTR but including 3'UTR
sub transcript_seq_spliced_without_5prime_utr {
   my ($tr) = @_;
   my $sequence = substr $tr->spliced_seq("soft_mask"), $tr->cdna_coding_start()-1, length $tr->spliced_seq("soft_mask");
   return $sequence;
}

#get transcript coding region
sub transcript_coding_seq {
   my ($tr) = @_;
   my $sequence = $tr->spliced_seq("soft_mask");
   $sequence =~ s/[a-z]//g; 
   return $sequence;
}

# slice sequence and insert variant of interest 
sub insert_variant_into_exon {
    my ($sequence, $ref_variant_seq, $variant, $exon_start, $exon_end, $var_start, $var_end) = @_;

    return "" if not defined($var_end);

    $variant =~ s/-//;
    my $before_variant = substr($sequence,0,$var_start - $exon_start);
    my $after_variant = substr($sequence, $var_end - $exon_start + 1, $exon_end - $var_end);
    my $variant_in_ref_seq = substr($sequence,$var_start - $exon_start,$var_end - $exon_start + 1 - $var_start + $exon_start); 
   
    return $before_variant . $variant . $after_variant;
}

#split seq in tri-grams denoting a list of codons
sub codon_ize_sequence {
    my ($sequence) = @_;
    $sequence =~ s/(\w{3})/$1 /g;
    return split(/ /,$sequence);
}

#determine whether condon is a stop-codon
sub is_stop_codon {
    my ($codon) = @_;
    $codon = uc $codon;
    return ($codon =~ /TAG/) || ($codon =~ /TAA/) || ($codon =~ /TGA/);
}

#pad the shorter of two sequences with additional charcters 
sub pad_sequence_pair {
    my ($seq_a, $seq_b) = @_;   
    my $diff = (length $seq_a) - (length $seq_b);
    $seq_a .= "X" x ($diff * -1) if $diff * -1;
    $seq_b .= "X" x ($diff) if $diff;
    return ($seq_a, $seq_b)
}

sub compare_codon_lists {
    my @codons_w_variant = @{$_[0]};
    my @codons_wo_variant = @{$_[1]};

    my $wt_sequence = "";
    my $added_sequence = "";
    my $stop_codon_cnt = 0;
    my $addition_started = 0;

    my $iter = each_array(@codons_wo_variant, @codons_w_variant);

    while( my($codon_a, $codon_b) = $iter->()) {

        # stop codon before variant is seen
        if(!$addition_started && is_stop_codon($codon_b)) {
            $stop_codon_cnt += 1;
            $termination_code = 2;
            last;
        }

        # if there is no codon change compared to reference sequence 
        if(!$addition_started && ($codon_a eq $codon_b)) {
            $wt_sequence .= $codon_a;
	    next;
        } else {
            $addition_started = 1;
        }

        # stop codon after variant
        if(is_stop_codon($codon_b)) {
            $stop_codon_cnt += 1;
            $termination_code = 3;
            last;
        } else {
            $added_sequence .= $codon_b;
        }
    }

    # if no stop codon was seen at all
    if($stop_codon_cnt < 1) {
        $termination_code = 4;
        return "","";
    }

    # final viable captured sequence
    return ($wt_sequence, $added_sequence);
}


sub process_transcript {
    my ( $transcript_variation_allele) = @_;
    my $transcript = $transcript_variation_allele->transcript();
    my $seq = transcript_seq_spliced_without_5prime_utr($transcript);
    my $ref_variant_seq = (split /\//, $transcript_variation_allele->allele_string)[0];
    $ref_variant_seq =~ s/-//;
    my $variant = $transcript_variation_allele->variation_feature_seq;

    my $condon_seq_with_variant = insert_variant_into_exon
    (
        $seq,
        $ref_variant_seq,
        $variant,
        $transcript->cdna_coding_start(),
        length $transcript->spliced_seq("soft_mask"),
        $transcript_variation_allele->transcript_variation->cdna_start,
        $transcript_variation_allele->transcript_variation->cdna_end
    );

    my $condon_seq_wo_variant = uc $seq;
    my ($condon_seq_with_variant_pad,$condon_seq_wo_variant_pad) = pad_sequence_pair($condon_seq_with_variant,$condon_seq_wo_variant);
    my @codons_w_var = codon_ize_sequence($condon_seq_with_variant_pad);
    my @codons_wo_var = codon_ize_sequence($condon_seq_wo_variant_pad);

    #return translate_seq_string( compare_codon_lists(\@codons_w_var, \@codons_wo_var));
    return compare_codon_lists(\@codons_w_var, \@codons_wo_var);
}

1;
