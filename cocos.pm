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
use Config::Simple;
use List::MoreUtils qw( pairwise );
use List::MoreUtils qw( each_array );

my $output_type = get_config('output_type');
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

sub get_config {
    my ($param) = shift;
    my $cfg = new Config::Simple('cocos.cfg');
    my $param_value = $cfg->param($param);
    return $param_value;
}

sub run {
    my ($self, $transcript_variation_allele, $vep_line) = @_;

    my $variant_cdna_pos = $vep_line->{'cDNA_position'};
    my $variant_cds_pos = $vep_line->{'CDS_position'};
    my $transcript_id = $vep_line->{'Feature'};
    my $uploaded_variation = $vep_line->{'Uploaded_variation'};
    $uploaded_variation =~ s/\//_/g;

    my $result = '';
    my $trans_len_pct_retained = 0;

    if($variant_cdna_pos !~ /^-$/ && $variant_cds_pos !~ /^-$/) {
        $variant_cdna_pos = (split /-/, $variant_cdna_pos)[0]; 

        $result = process_transcript( $transcript_variation_allele, $variant_cdna_pos);

        # get chr and pos for current variant
        my $chr = $transcript_variation_allele->transcript->feature_Slice->seq_region_name;
        my $pos = ($transcript_variation_allele->transcript->get_TranscriptMapper->cdna2genomic($variant_cdna_pos,$variant_cdna_pos))[0]->start;
	my $chr_pos = $chr."_".$pos;
        my $allele_str = $transcript_variation_allele->allele_string;

        if($result ne '') {
            write_to_file($output_type, $chr, $pos, get_config('output_path'), $transcript_id, $uploaded_variation, $allele_str, $result);
        }

	my $trans = $transcript_variation_allele->transcript;
        $trans_len_pct_retained = ($variant_cdna_pos - $trans->cdna_coding_start)/($trans->cdna_coding_end - $trans->cdna_coding_start);	
	    
        return {
            TRANS_LEN => ($trans_len_pct_retained > 0 && length $result > 0 ? $trans->length : "N/A"),
            TRANS_RET_PCT => ($trans_len_pct_retained > 0 && length $result > 0 ? substr($trans_len_pct_retained,0,7) : "N/A"),
            ALT_AA_SEQ => (length $result > 0 ? length $result : 'N/A'),
	    TERM_TYPE => $termination_code,
	    FILE_NAME => $chr_pos
        };
    }

    return {};
}

sub write_to_file {
    my ($output_type, $chr, $pos, $file_path, $trans_id, $uploaded_var, $allele_str, $seq) = @_;

    my $content = '';
    my $filename = '';
    my $seq_len = length $seq;

    if($output_type eq 'fasta') {
	# default fasta length of 80 characters per line
        my $len = 80;

    	$content .= ">$chr $pos $trans_id $seq_len $uploaded_var $allele_str\n";
    	while (my $chunk = substr($seq, 0, $len, "")) {
            $content .= "$chunk\n";
    	}
 
        $filename = $file_path.'/'.$chr.'_'.$pos.".fasta";
    }

    if($output_type eq 'tsv') {
        $filename = $file_path.'/cocos_results.tsv'; 
	$content .= $chr."\t".$pos."\t".$trans_id."\t".$uploaded_var."\t".$allele_str."\t".$seq_len."\t".$seq."\n";
    }

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


# iterate through all exons and find exon where variant is located in a transcript
sub process_transcript {
    my ($transcript_variation_allele, $variant_cdna_pos) = @_;
    return process_sequence( $transcript_variation_allele); 
}

sub transcript_seq_spliced_without_5prime_utr {

   my ($tr) = @_;
   my $sequence = substr $tr->spliced_seq("soft_mask"), $tr->cdna_coding_start()-1, length $tr->spliced_seq("soft_mask"); 
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

sub codon_ize_sequence {
    my ($sequence) = @_;
    $sequence =~ s/(\w{3})/$1 /g;
    return split(/ /,$sequence);
}


sub is_stop_codon {
    my ($codon) = @_;
    $codon = uc $codon;
    return ($codon =~ /TAG/) || ($codon =~ /TAA/) || ($codon =~ /TGA/);
}

sub pad_sequence_pair {
    my ($seq_a, $seq_b) = @_;
    
    my $diff = (length $seq_a) - (length $seq_b);
    $seq_a .= "X" x ($diff*-1) if $diff*-1;
    $seq_b .= "X" x ($diff) if $diff;

    return ($seq_a, $seq_b)
}

sub compare_codon_lists {
    my @codons_w_variant = @{$_[0]};
    my @codons_wo_variant = @{$_[1]};

    my $added_sequence = "";
    my $stop_codon_cnt = 0;
    my $addition_started = 0;

    my $iter = each_array(@codons_wo_variant, @codons_w_variant);

    while ( my ($codon_a, $codon_b) = $iter->()) {

        # stop codon before variant is seen
        if(!$addition_started && is_stop_codon($codon_b)) {
            $stop_codon_cnt += 1;
            $termination_code = 2;
            last;
        }

        # if there is no codon change compared to reference sequence 
        if(!$addition_started && ($codon_a eq $codon_b)) {
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

    # if no stop codon was seen 
    if($stop_codon_cnt < 1) {
        $termination_code = 4;
        return "";
    }

    # final viable captured sequence
    return $added_sequence;
}


sub process_sequence {
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


    my $condon_seq_wo_variant = uc $seq;#= initialize_exon_sequence($exon);

    my ($condon_seq_with_variant_pad,$condon_seq_wo_variant_pad) = pad_sequence_pair($condon_seq_with_variant,$condon_seq_wo_variant);

    my @codons_w_var = codon_ize_sequence($condon_seq_with_variant_pad);
    my @codons_wo_var = codon_ize_sequence($condon_seq_wo_variant_pad);

    return translate_seq_string( compare_codon_lists(\@codons_w_var, \@codons_wo_var));
}

1;
