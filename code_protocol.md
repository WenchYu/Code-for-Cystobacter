
# Figure 1 genomic features
## Figure 1A phelogentic tree
### 16S 建树
working directory: ws2 /media/weibin/sde2/cystobacter/prokka_results
1. 提取16S
```bash
find . -type f -name "*.fna" > ../mono.txt 
prokka_prediction
./do_prokka.sh
```
2. 从 fnn DNA注释文件中挑出最长的 16S ribosomal RNA 汇总成一个 fasta
3. 建树
    `kalign -i cystobacter_16S.fna -f fasta > 16s-align.fasta` #记得除去 *align* 文件开头的注释
    `fasttree -nt 16s-align.fasta > cystobacter_tree.txt `#获取树文件
### UBCG 建树
working directory: hehe/
1. 切换 java
export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.8.0_291.jdk/Contents/Home
export PATH=$JAVA_HOME/bin:$PATH
java -version
2. 提取 bcgs 
run_ubcg.sh
3. align
java -jar UBCG.jar align -bcg_dir bcg -prefix cystobacter
### checkm 
查看 completeness 和 contamination
checkm lineage_wf -t 20 -x fna --nt --tab_table -f cystobacter.txt fasta/ results-cystobacter/
lineage_wf：运行 CheckM 的完整工作流（workflow）。
-t 20：使用 20 线程进行计算。
-x fna：输入文件的扩展名为 .fna。
--nt 生成基因序列
--tab_table：输出格式为制表符分隔的表格。
-f cystobacter.txt：将结果写入 dsm43850.txt。
fasta/：输入基因组目录（你的 fna 文件应该放在这里）。
results-cystobacter/：存放 CheckM 计算结果的目录。
### seqkit 查看GC content
seqkit stats --all fasta/*.fna

## Figure 1B Similarities of genome

[JSspeciesWS](https://jspecies.ribohost.com/jspeciesws)
The results can be analysed with the sharing code: 37C2DE80CC9F94A142F9

## Fig 1C BGC classification
### bigscape
`
nohup python /media/weibin/5294B75294B736F7/BiG-SCAPE-1.1.7/bigscape.py -c 36 --pfam_dir /media/weibin/5294B75294B736F7/Pfam -i /media/weibin/W16/ywc/gbk -o /media/weibin/W16/uwc/bigscape --mibig --mode global --mix --include_gbk_str "*" --include_singletons --cutoffs 0.5 0.7 &
`

-i gbk/ 其中是不包括汇总gbk的所有antismash预测的.gbk文件

###  
Network_Annotations_Full.tsv中包含对各genome的汇总

# Figure 2
## MCS similarity network
`MCS-SN.py` -> `MCS-SN_Sim0.7Topk10.cys`
activity summary [CMNPD活性分类](https://docs.cmnpd.org/frequently-asked-questions/data-questions#what-do-the-abbreviations-in-biological-activity-data-mean)

## BGC sequence similarity network

Archangium violaceum strain Cb vi35 (formerly known as Cystobacter violaceus strain Cb vi35)
Archangium violaceum strain Cb vi105 (formerly known as Cystobacter violaceus strain Cb vi105)
Archangium violaceum strain Cb vi76 (formerly known as Cystobacter violaceus strain Cb vi76)

BGC0001413 --> Cystobacter sp. Cbv34 Cystobactamides
BGC0000335 --> Cystobacter fuscus DSM 2262 Cystomanamides
BGC0001344  --> Cystobacter sp. SBCb004 Tubulysins

# Figure 3 antimicobial peptide
PyCharm -> settings -> Project -> project structure
添加 astool 和 script 为 sources
add path也没什么用
copy 到新的脚本里倒是可以用

## extract ripps
`ex_ripps.py` -> 'ripps_cystobacter.tsv'
直接用总的 `_genomic.gbk` 提取， 看不到具体的 region几
1. 汇总所有 region_xxx.gbk文件

## ripps sturctural similarity network
揭示短肽之间的关系： 若关注功能：优先理化性质或k-mer方法； 若关注进化关系：选择序列比对或语言模型。
`tsv2fast.py` -> 'ripps_cystobacter.fasta'
antismash 预测的时候, 相同位置保留最高值， 不同位置全部保留
`KmerSSN.py` -> 'CorePeptide_K2Sim0.8Top5.graphml'

## antimicrobial prediction
[AMP Scanner](https://www.dveltri.com/ascan/v2/index.html)
[camp](http://www.camp3.bicnirrh.res.in/index.php)

>GCA_044360605.1——C. fuscus Cbf8——CP043492.1.region037.gbk_1 28.6% against NZ_FNGC01000009.1
TRSTARCFCTSSSTCTVSGACPHRASRPGKVWW
>GCA_044360605.1——C. fuscus Cbf8——CP043492.1.region037.gbk_3	
SCGDSRVCVGDCCADGV
>GCA_044360685.1——C. fuscus Cbf10——CP043491.1.region007.gbk_2 31.6% against NZ_CP067395.1
LSKKVWCAPGGCDGKS
>GCA_002305875.1——C. fuscus DSM 52655——CP022098.1.region007.gbk_6 25.8% against NZ_CP069338.1
LSPDRWRKRRHLAARSRRCWPCGW 

>GCA_000335475.2——C. fuscus DSM 2262——ANAH02000066.1.region002.gbk_4 23.3% against NZ_JAHFVJ010000001.1
TTARSLCTCALPPQVFSSPPGRIARILRHPLRPTASKLPIPSH
ANAH02000066.1.region002.gbk_1 50.0% against NZ_KB913037.1
AGVNASCGWSSCNRTN


>GCA_002305875.1——C. fuscus DSM 52655——CP022098.1.region007.gbk_2_1.7 12.1% against NZ_FOFV01000025.1
HPRTEHFKEDIIQAFYDGIKHKPRTTFGNVKADVIADKEPLFIRGNFCRVIRESAWRG


>GCA_000335475.2——C. fuscus DSM 2262——ANAH02000014.1.region001.gbk_1 91.3% against NZ_CP022098.1/ 14.3% against NC_019757.1
RCVRCHAREHMRTSHRVAGRRRS


## Seqlogo
可视化 Core peptide

## Weblogo3
Kalign
>CP043492.1.region037.gbk_1
TRSTARCFCTSSSTCTVSGACPHRASRPGKVWW
>CP043491.1.region007.gbk_8
ARHTASHPVSSSAAPSGRCHVGTMAVAAKRSVLMP
>CP022098.1.region007.gbk_1
YKYNRRQKQWMSWWGAGSGGVHFRQGMDLECADEHRCGCFVGMLWPAWCV
>ANAH02000066.1.region002.gbk_4
TTARSLCTCALPPQVFSSPPGRIARILRHPLRPTASKLPIPSH
>CP022098.1.region007.gbk_5
VVEAGSFSAAARALVMNSGALRGLPTPCTQANARSNKTKPFRKDDRPAGI
>CP043492.1.region010.gbk_1
LGLWATAAWRSTPCSLAPGVGVAARARPTSTELSSSSVTSWKNSSG
>CP022098.1.region039.gbk_2
WKDPDPVQHSCHPSPPRTHERPILRLFQDVFQGVSRKPFTLRTSSQALTGL
>CP022098.1.region007.gbk_3
PAGIVELSDEALDSLLSGGKQVSSCCWESC
>CP022098.1.region007.gbk_2
HPRTEHFKEDIIQAFYDGIKHKPRTTFGNVKADVIADKEPLFIRGNFCRVIRESAWRG
>CP043492.1.region010.gbk_6
ASGLGTPYRLPAAGPSTWRPAPASRPRVRWRARGRDRVRTSASSRSPCAA
>ANAH02000066.1.region002.gbk_3
PSTLGGRQPHAPSAPALSLHRSSLRRQVASLEFFDTLSGLPQANCQFLLIERLRRGA
>CP043492.1.region037.gbk_2
QALSTSCGPVLTLSSGCYIKLDLDEPYVEALTYAAPSSAVIRRSGGPPGSPLP


## Ripps BGCs SSN
GCA_036387055.1 C. sp. SMAG_U3318 和 GCA_913775605.1 C. sp. SP168 组装质量不好考虑去掉
直接从总的里面挑出