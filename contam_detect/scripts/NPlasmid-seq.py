# -*- coding: utf-8 -*-
# @Time    : 2022/12/07 9:50
# @Author  : zzk
# @File    : main.py
# Description:
import argparse
import os

import FromSamExtract
import SimpleSam
import WaringCluster
import consensus
import Visualize
import MakeNewRef
import Report
import glob
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')
    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('-d', '--data',
                            help="the Sam file 、 Reference file and DemultiplexInfo file(the same name)", required=True)
    all_parser.add_argument('-b', '--baseSequence', default=None, help="the unique 20 base")
    all_parser.add_argument('-n', '--num', default=None, help="Cluster the first few classes")

    visualize_parser = subparsers.add_parser('visualize',
                                             help='Visual sam generation for fastq files')
    visualize_parser.add_argument('-d', '--data', help="the Reference file and DemultiplexInfo file(the same name)",
                                  required=True)

    extract_parser = subparsers.add_parser('extract',
                                           help='from fastq to extract waring "sgRNA 20 nt" to a new fastq')
    extract_parser.add_argument('-d', '--data', help="the Reference file and DemultiplexInfo file(the same name)",
                                required=True)
    extract_parser.add_argument('-b', '--baseSequence', help="the unique 20 base", required=True)

    cluster_parser = subparsers.add_parser('cluster',
                                           help='Put the waring-"sgRNA 20 nt" clustering')
    cluster_parser.add_argument('-w', '--waring_fq', help="the waring data file", required=True)
    cluster_parser.add_argument('-n', '--num', help="Cluster the first few classes", required=True)

    return parser.parse_args()


'''
2022-12-20 现在整体是sgRNA污染，从原始fastq和原始fasta到提异常数据的共识序列。
'''
import docx
if __name__ == '__main__':
    args = parse_args()
    if not os.path.exists("./consen"):
        os.mkdir("./consen")
    if args.command == 'all':
        logger.info("start")
        path = f'{args.data}'
        seq_20nt = args.baseSequence
        if seq_20nt:
            top_n = int(args.num)
            fasta_file = glob.glob(f"{path}/*.fasta")
            for file_fa in fasta_file:  # 三个原始文件，file_fa,file_fq,file_sam
                file_fa, lef_seq = MakeNewRef.run(file_fa=file_fa, basic_20=seq_20nt)
                Visualize.run(file=file_fa)
                file_fq, file_sam = f"{file_fa[:-6]}.fastq", SimpleSam.run(f"{file_fa[:-6]}.sam")
                logger.info(f"Start to extract warning data({file_fq})")
                waring_fq, allData = FromSamExtract.main(file_fa=file_fa, file_fq=file_fq, file_sam=file_sam,
                                                         seq_20nt=seq_20nt)
                WaringCluster.run(waring_fq=waring_fq, top_n=top_n)
                dic_basic20 = consensus.run(left_seq=lef_seq,file_fa=file_fa)  # 共识序列保存在./racon/consensus_racon.fa
                # print(dic_basic20)
                Report.run(allData=allData, dic_basic20=dic_basic20, fa_file=f'{file_fa[:-6]}-ori.fasta',
                           seq_20nt=seq_20nt)

        else:
            '''
            2022.12.21普通污染判断
            '''
            fasta_file = glob.glob(f"{path}/*.fasta")
            allData = {}
            cluster_fastq = {}  # fasta名做key，[]做value
            for file_fa in fasta_file:  # 三个原始文件，file_fa,file_fq,file_sam
                cluster_fastq[f'{file_fa}'] = []
                Visualize.run(file=file_fa)
                file_fq, file_sam = f"{file_fa[:-6]}.fastq", SimpleSam.run(f"{file_fa[:-6]}.sam")
                logger.info(f"Start to extract warning data({file_fq})")
                waring_fq, waring_blast_result, allData[f'{file_fa}'] = FromSamExtract.main(file_fa=file_fa,
                                                                                            file_fq=file_fq,
                                                                                            file_sam=file_sam,
                                                                                            seq_20nt=None,
                                                                                            refervalue=0.9,
                                                                                            readvalue=0.9)  # fq、fa、sam
                cluster_fastq[f'{file_fa}'] = WaringCluster.commonRun(waring_fq=waring_fq,
                                                                      waring_blast_result=waring_blast_result,
                                                                      top_n=int(args.num))
            dic_basic20 = consensus.consensus()  # 共识序列保存在./racon/consensus_racon.fa
            for file_fa in fasta_file:
                doc = docx.Document()
                doc_file = Report.run(doc=doc, allData=allData[f'{file_fa}'], fa_file=file_fa, seq_20nt=None,
                                      dic_basic20=None,
                                      top_n_file=cluster_fastq[f'{file_fa}'])
                '''
                Visualize.run(file=file_fa)
                file_fq, file_sam = f"{file_fa[:-6]}.fastq", SimpleSam.run(f"{file_fa[:-6]}.sam")
                logger.info(f"Start to extract warning data({file_fq})")
                waring_fq, waring_blast_result, allData = FromSamExtract.main(file_fa=file_fa, file_fq=file_fq,
                                                                              file_sam=file_sam, seq_20nt=None,
                                                                              refervalue=0.8,
                                                                              readvalue=0.8)  # fq、fa、sam
                WaringCluster.commonRun(waring_fq=waring_fq, waring_blast_result=waring_blast_result,
                                        top_n=int(args.num))
                dic_basic20 = consensus.consensus()  # 共识序列保存在./racon/consensus_racon.fa
                Report.run(allData=allData, fa_file=file_fa,seq_20nt=None,dic_basic20=None)

                '''


    elif args.command == 'visualize':
        path = f'{args.data}'
        fasta_file = glob.glob(f"{path}/*.fasta")
        for file_fa in fasta_file:  # 三个原始文件，file_fa,file_fq,file_sam
            Visualize.run(file_fa)

    elif args.command == 'extract':
        path = f'{args.data}'
        seq_20nt = args.baseSequence
        fasta_file = glob.glob(f"{path}/*.fasta")
        for file_fa in fasta_file:  # 三个原始文件，file_fa,file_fq,file_sam
            file_fq, file_sam = f"{file_fa[:-6]}.fastq", SimpleSam.run(f"{file_fa[:-6]}.sam")
            logger.info(f"Start to extract warning data({file_fq})")
            waring_fq = FromSamExtract.main(file_fa=file_fa, file_fq=file_fq, file_sam=file_sam, seq_20nt=seq_20nt)

    elif args.command == 'cluster':
        waring_fq = args.waring_fq
        logger.info(f"Start clustering warning data({waring_fq})")
        top_n = int(args.num)
        WaringCluster.run(waring_fq=waring_fq, top_n=top_n)
        # consensus.run()
