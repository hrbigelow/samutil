            //tally running depth measurement
            guide_depth += prev_feature.guide_jumps;
            correct_depth += prev_feature.correct_jumps;
            error_depth += prev_feature.error_jumps;

            assert(guide_depth > 0);
            assert(correct_depth >= 0);
            assert(error_depth >= 0);

            //could define regions based on thresholds, or print the coverage
            //plot.  
            
            //print it
            fprintf(output_depth_fh, "%s\t%Zu\t%Zu\t%Zu\t%Zu\n",
                    prev_contig, prev_feature_bound + 1, 
                    guide_depth, correct_depth, error_depth);



    // size_t contig_bound;

    // feature const* current_feature;
    // for (OFFSET_ITER ci = contig_offsets.begin(); 
    //      ci != contig_offsets.end(); ++ci)
    // {
    //     char const* contig_name = (*ci).first;
    //     size_t contig_offset = (*ci).second;
    //     size_t contig_length = contig_lengths[(*ci).first];

    //     for (current_feature = eval_buffer + contig_offset, contig_bound = 0;
    //          contig_bound != contig_length; ++contig_bound, ++current_feature)
    //     {
    //         guide_depth += current_feature->guide_jumps;
    //         correct_depth += current_feature->correct_jumps;
    //         error_depth += current_feature->error_jumps;

    //         assert(guide_depth >= 0);
    //         assert(correct_depth >= 0);
    //         assert(error_depth >= 0);

    //         if (guide_depth > 0 || correct_depth > 0 || error_depth > 0)
    //         {
    //             fprintf(output_depth_fh, "%s\t%Zu\t%Zu\t%Zu\t%Zu\n",
    //                     contig_name, contig_bound + 1, 
    //                     guide_depth, correct_depth, error_depth);
    //         }
            
    //     }
        
    // }
