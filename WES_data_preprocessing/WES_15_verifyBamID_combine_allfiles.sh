#!/bin/bash

head -1 /WES_output_combine/verifyBamID_files/R0218_L005_verifyBamID.selfSM > /WES_output_combine/verifyBamID_files/All_Samples_VerifyBamID.selfSM
for i in /WES_output_combine/verifyBamID_files/*.selfSM
do
tail -1 $i >> /WES_output_combine/verifyBamID_files/All_Samples_VerifyBamID.selfSM
done

exit 0
