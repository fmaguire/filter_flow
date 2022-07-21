#!/usr/bin/env python

from logging import exception
import re
import os
import argparse
import shutil
import pysam
import gzip
import datetime
import json

#structure for sotring snps
class CARD:
    class ARO:
        def __init__(self, type, data, index = None) :
            data = self.Parse(type, data, index)
            self.aro = data["ARO"]
            self.name = data["Name"]
            self.family = data["Family"]
            self.species = data["Species"]
            self.proteinCluster = ""
            self.dnaCluster = ""
            if ("DNA" in data.keys()):
                self.dnaAccession = data["DNAAccession"]
                self.DNA = data["DNA"]
                self.proteinAccession = None
                self.AA = None
            if ("AA" in data.keys()):
                self.proteinAccession = data["ProteinAccession"]
                self.AA = data["AA"]
                self.dnaAccession = None
                self.DNA = None         
    
        def __str__(self):
            return str(self.aro) + ":" + self.name
        
        def Parse(self, type, input, index = None):
            #nucl: >gb|GQ343019|+|132-1023|ARO:3002999|CblA-1 [mixed culture bacterium AX_gF3SD01_15] 
            if (type == "DNA"):
                headers = "gb\tDNAAccession\tStrand\tPosition\tARO\tName\tDNA"
            elif (type == "AA"):
                headers = "gb\tProteinAccession\tARO\tName\tAA"
            else:
                raise exception("CARD DB type is neither DNA or AA")
            metadata = {}
            input = input.split("|")
            headers = headers.split("\t")
            col = 0
            for header in headers:
                metadata[header] = input[col].strip()
                col = col + 1
            metadata["Species"] = metadata["Name"].split("[")[1].replace("]","").strip()
            metadata["Name"] = metadata["Name"].split("[")[0].strip()
            metadata["ARO"] = metadata["ARO"].replace("ARO:","")
            metadata["Family"] = index[metadata["ARO"]]["AMR Gene Family"]
            return metadata
        
        def GetARO(self):
            return self.aro


    def __init__(self, dnaPath, aaPath, indexPath) -> None:
        self.index = self.ParseIndex(indexPath)
        self.CARD_AA = self.LoadProteinData(aaPath, self.index)
        self.CARD_DNA = self.LoadDNAData(dnaPath, self.index)
        self.CARD_All = self.CombineData(self.CARD_DNA, self.CARD_AA)
        self.clusterDict = None
        self.accessionDict = None
    
    def ParseIndex(self, indexPath): #discards metamodels
        headers = "ARO Accession\tCVTERM ID\tModel Sequence ID\tModel ID\tModel Name\tARO Name\tProtein Accession\tDNA Accession\tAMR Gene Family\tDrug Class\tResistance Mechanism"        
        index = {}
        headerList = headers.split("\t")


        with open(indexPath, 'r') as fh:
            for line in fh:
                if (line != headers):
                    metadata = {}
                    input = line.strip().split("\t")
                    col = 0
                    for header in headerList:
                        metadata[header] = input[col].strip()
                        col = col + 1
        
                    metadata["ARO Accession"] = metadata["ARO Accession"].replace("ARO:","")
                    index[metadata["ARO Accession"]] = metadata

        return index

    
    def LoadDNAData(self, path, index):
        card = [line.strip() for line in open(path, 'r')]
        aro = {}
        for i in range(0, len(card), 2):
            data = card[i] + "|" + card[i+1]
            aroObj = self.ARO("DNA", data, index)
            aro[aroObj.GetARO()] = aroObj
        return aro
    
    def LoadProteinData(self, path, index):
        card = [line.strip() for line in open(path, 'r')]
        aro = {}
        for i in range(0, len(card), 2):
            data = card[i] + "|" + card[i+1]
            aroObj = self.ARO("AA", data, index)
            aro[aroObj.GetARO()] = aroObj
        return aro
    
    def CombineData(self, dnaData, proteinData):
        keys = set(list(dnaData.keys()) + list(proteinData.keys()))
        aro = {}
        for key in keys:
            if (key in dnaData.keys()):
                aro[key] = dnaData[key]
                if (key in proteinData.keys()):
                    aro[key].proteinAccession = proteinData[key].proteinAccession
                    aro[key].AA = proteinData[key].AA
            elif (key in dnaData.keys()):
                aro[key] = dnaData[key]
                if (key in dnaData.keys):
                    if aro[key].dnaAccession == None:
                        aro[key].dnaAccession = dnaData[key].dnaAccession
                        aro[key].DNA = dnaData[key].DNA
        return aro


    def GetCARDDNA(self):
        return self.CARD_DNA
    
    def GetCARDAA(self):
        return self.CARD_AA
    
    def GetCARDALL(self):
        return self.CARD_All         

    def GetClusterDict(self):
        return self.clusterDict   
    
    def GetAccessionDict(self):
        return self.accessionDict   

    def AssignHMMClusters(self, dnaClusterPath, proteinClusterPath):
        card = self.CARD_All

        for i in [file for file in os.listdir(dnaClusterPath) if file.endswith(".msa_clean") ]:
            clusterName = (i.split(".")[0])
            aroInCluster = set([line.split("ARO:")[1].split("|")[0] for line in open(dnaClusterPath + "/" + i,'r') if line.find("ARO:") != -1])
            for aro in aroInCluster:
                if card[aro].dnaCluster != "":
                    raise exception("cluster already assigned in nucl")
                card[aro].dnaCluster = clusterName
        for i in [file for file in os.listdir(proteinClusterPath) if file.endswith(".faa") ]:
            clusterName = (i)
            aroInCluster = set([line.split("ARO:")[1].split("|")[0] for line in open(proteinClusterPath + "/" + i,'r') if line.find("ARO:") != -1])
            for aro in aroInCluster:
                if card[aro].proteinCluster != "":
                    raise exception("cluster already assigned in prot")
                card[aro].proteinCluster = clusterName

        clusterDict = {}
        accessionDict = {}
        for aro in card.values():
            if (aro.proteinAccession != ""):
                accessionDict[aro.proteinAccession] = aro.aro
            if (aro.dnaAccession != ""):
                accessionDict[aro.dnaAccession] = aro.aro
            if (aro.proteinCluster != ""):
                clusterDict[aro.proteinCluster] = aro.aro
            if (aro.dnaCluster != ""):
                clusterDict[aro.dnaCluster] = aro.aro
        self.clusterDict = clusterDict
        self.accessionDict = accessionDict

class Results:
    def __init__(self, read, aro, name, tool, headers, metadata):
        self.read = read
        self.aro = aro
        self.name = name
        self.tool = tool
        self.headers = headers
        self.metadata = metadata
        self.truth = None
        self.aroFamily = None
        self.truthAROFamily = None

    def __str__(self) -> str:
        return ("{}\t{}\t{}\t{}\t{}\t{}".format(self.read, self.aro, self.name, self.aroFamily, self.truth, self.truthAROFamily))
    
    def __eq__(self, __o: object) -> bool:
        if (isinstance(__o, Results)):
            return (self.read == __o.read and 
                self.aro == __o.aro and 
                self.tool == __o.tool)
        return False
    
    def __hash__(self) -> int:
        return hash(self.read + str(self.aro) + self.tool)

    def AppendTruthARO(self, resultObj, truthObjDict):
        resultObjReadName = resultObj.read.split("/")[0]
        if (resultObjReadName in truthObjDict.keys()):
            resultObj.truth = truthObjDict[resultObjReadName].aro
        else:
            resultObj.truth = "FP"
        return resultObj
    
    def AppendAROFamily(self, resultObj, cardObj):
        resultObj.aroFamily = ";".join([cardObj.GetCARDALL()[a].family for a in resultObj.aro.split(";")])
        if (resultObj.truth == "FP"):
            return resultObj
        resultObj.truthAROFamily = ";".join([cardObj.GetCARDALL()[a].family for a in resultObj.truth.split(";")])
        return resultObj

    def Parse(self):
        raise exception("base result class have no defined function for Parse, use a child")
    
    def Filter(self,resultList):
        return resultList
    
    def CleanReadName(self, read):
        if (read.find("/") != -1):
            if (read.split("/")[1].find("_") != -1):
                #this is likely a protein
                cleanRead = read.split("/")[0] + "/" + read.split("/")[1].split("_")[0]
            else:
            #this is nucl, 
                cleanRead = read
        else:
            cleanRead = read
        return cleanRead
    
    def GetData(self):
        #this should return a list of cleaned data of format: Read > Label
        return self.metadata
    
    def GetHeaders(self):
        return self.headers    

#data class for outfmt6 results. i.e. diamond/blast
class Fmt6Results(Results):
    def __init__(self, tool, rawData, headers= None,metadata = None, card=None):
        if (metadata == None):
            headers = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore"
            data = self.Parse(headers, rawData)
        else:
            headers = headers
            data = metadata
        self.read = data["qseqid"]
        self.aro = data["sseqid"].split("|")[4].split(":")[1]
        self.name = data["sseqid"].split("|")[5]
        self.pident = data["pident"]
        self.length = data["length"]
        self.mismatch = data["mismatch"]
        self.gapopen = data["gapopen"]
        self.qstart = data["qstart"]
        self.qend = data["qend"]
        self.sstart = data["sstart"]
        self.send = data["send"]
        self.evalue = data["evalue"]
        self.bitscore = data["bitscore"]
        self.metadata = metadata

        super().__init__(self.read, self.aro, self.name, tool, headers, data)
    
    def Parse(self, headers, rawData):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        #returns a dictionary of columnName:data
        metadata = {}
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        col = 0
        for header in headers:
            metadata[header] = data[col].strip()
            col = col + 1
        
        return metadata

    def Filter(self, resultList, newReadObj, evalue = 1e-10, bitscore = 100):
        newReadName = self.CleanReadName(newReadObj.read)
        for i in range(0, len(resultList)):
            if float(resultList[i].evalue) > evalue and float(resultList[i].bitscore) < bitscore:
                #remove the entry cuz its below the cutoff
                del resultList[i]
            else:
                if (len(resultList) == 0):
                    if float(resultList[i].evalue) < evalue or float(resultList[i].bitscore) > bitscore:
                        resultList.append(newReadObj)
                else:
                    uniqueRead = self.CleanReadName(resultList[i].read)
                    if (uniqueRead == newReadName):
                        #a new hit with the same read name. i.e. read123/1 
                        if float(newReadObj.evalue) < float(resultList[i].evalue) or float(newReadObj.bitscore) > float(resultList[i].bitscore):
                            resultList[i] = newReadObj
                    else:
                        print("test")
                        resultList.append(newReadObj)
        return resultList

class AAFmt6Result(Fmt6Results):
    def __init__(self, tool, rawData, headers= None,metadata = None, card=None):
        if (metadata == None):
            headers = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore"
            data = self.Parse(headers, rawData)
        else:
            headers = headers
            data = metadata

        super().__init__(tool, rawData,headers, data)
    
    def Parse(self, headers, rawData):
        #there are a few extra | in the nucleotide fasta but not protein fasta. just add them then call super.
        rawData = rawData.replace("|ARO:", "|||ARO:")
        return super().Parse(headers, rawData)
    
    def Filter(self, resultList, newReadObj):
        return super().Filter(resultList, newReadObj, 1e-10, 100)

#data class for mmseq2's extended outfmt6 results
class MMSeqFmt6Result(Fmt6Results):
    def __init__(self, tool, rawData, card):
        headers = "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qlen	slen"
        data = self.Parse(headers, rawData, card)
        self.qlen = data['qlen']
        self.slen = data['slen']
        super().__init__(tool, rawData, headers, data)
    
    def Parse(self, headers, rawData, card):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore, qlen, slen
        #returns a dictionary of columnName:data
        metadata = {}
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        headers = headers[:-2]
        metadata["qlen"] = data[-2].strip()
        metadata["slen"] = data[-1].strip()
        data = '\t'.join(data[:-2])

        metadata.update(super().Parse(headers, data))

        #format the seqids into the same format as the parent class
        #temp solution until i can format it properly
        #have to pass in a cluster>ARO mapping file path
        ARO = card.GetAccessionDict()[metadata["sseqid"]]
        data = card.GetCARDALL()[ARO]
        metadata["sseqid"] = "gb|"+ metadata["sseqid"] + "|||ARO:" + ARO + "|" + data.name

        return metadata

#data class for hmmsearch's tabular results
class HMMResults(Results):
    def __init__(self, tool, rawData, card=None):
        headers = "target name	accession	query name	accession	E-value	score	bias	E-value	score	bias	exp	reg	clu	ov	env	dom	rep	inc	description of target"
        data = self.Parse(headers, rawData, card)
        self.read = data["target name"]
        self.accession = data["accession"]
        self.queryname = data["query name"]
        self.aro = self.queryname.split("|")[0] 
        self.name = self.queryname.split("|")[1]
        self.evalue = data["E-value"]
        self.score = data["score"]
        self.bias = data["bias"]
        self.exp = data["exp"]
        self.reg= data["reg"]
        self.clu= data["clu"]
        self.ov= data["ov"]
        self.env= data["env"]
        self.dom= data["dom"]
        self.rep= data["rep"]
        self.inc = data["inc"]
        self.description = data["description of target"]

        super().__init__(self.read, self.aro, self.name, tool, headers, data)
    
    def Parse(self, headers, rawData, card):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: target name,accession  query name,accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #returns a dictionary of columnName:data
        metadata = {}
        rawData = re.sub(r"\s+", '\t', rawData)
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        col = 0
        for header in headers:
            metadata[header] = data[col].strip()
            col = col + 1

        ARO = card.GetClusterDict()[metadata["query name"]] 
        data = card.GetCARDALL()[ARO]
        metadata["query name"] = ""+ ARO + "|" + data.name

        return metadata
    
    def Filter(self, resultList, newReadObj, evalue = 1e-10, bitscore = 100):
        newReadName = self.CleanReadName(newReadObj.read)
        for i in range(0, len(resultList)):
            if float(resultList[i].evalue) > evalue and float(resultList[i].score) < bitscore:
                #remove the entry cuz its below the cutoff
                del resultList[i]
            else:
                if (len(resultList) == 0):
                    if float(resultList[i].evalue) < evalue or float(resultList[i].score) > bitscore:
                        resultList.append(newReadObj)
                else:
                    uniqueRead = self.CleanReadName(resultList[i].read)
                    if (uniqueRead == newReadName):
                        #a new hit with the same read name. i.e. read123/1 
                        if float(newReadObj.evalue) < float(resultList[i].evalue) or float(newReadObj.score) > float(resultList[i].score):
                            resultList[i] = newReadObj
                    else:
                        print("test")
                        resultList.append(newReadObj)
        return resultList

#data class for SAM/BAM results. i.e. bwa/bowtie2
class FmtSAMResults(Results):
    def __init__(self, tool, rawData, headers = None, metadata = None, card=None):
        if (metadata == None):
            headers = "ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual	Optionals"
            data = self.Parse(headers, rawData)
        else:
            headers = headers
            data = metadata

        self.read = data["ReadName"]
        self.flags = data["Flags"]
        self.refName = data["RefName"]
        if (self.refName != "*"):
            self.aro = data["RefName"].split("|")[4].split(":")[1]
            self.name = data["RefName"].split("|")[5]
        else:
            self.aro = "-1"
            self.name = "*"
        self.ForwardBase = data["ForwardBase"]
        self.MapQual = data["MapQual"]
        self.CIGAR = data["CIGAR"]
        self.MateRef = data["MateRef"]
        self.MateBase= data["MateBase"]
        self.FragLength= data["FragLength"]
        self.ReadSeq= data["ReadSeq"]
        self.ReadQual= data["ReadQual"]
        self.Optionals= data["Optionals"]

        self.metadata = metadata

        super().__init__(self.read, self.aro, self.name, tool, headers, data)
    
    def Parse(self, headers, rawData):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual	Optionals
        #returns a dictionary of columnName:data
        metadata = {}
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        col = 0
        for header in headers:
            if (col == 11):
                metadata[header] = "\t".join(data[col:])
                break
            metadata[header] = data[col].strip()
            col = col + 1
        return metadata

    def Filter(self, resultList, newReadObj, evalue = 1e-10, bitscore = 100):
        #resultList is a list of AAFmt6Result Objects
        result = [r for r in resultList if r.refName != "*"]
        for i in range(0, len(result)):
            if (result[i].read == newReadObj.read):
                if (newReadObj.refName != "*"):
                    result[i].aro = ";".join(set([result[i].aro, newReadObj.aro]))
                    result[i].name = ";".join(set([result[i].name, newReadObj.name]))
        return result

class GrootFmtSAMResults(FmtSAMResults):
    def __init__(self, tool, rawData, card=None):
        headers = "ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual	Optionals"
        data = self.Parse(headers, rawData)
        super().__init__(tool,  rawData, headers, data)
    
    def Parse(self, headers, rawData):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual Optionals
        #returns a dictionary of columnName:data
        #this is missing the optional field, so just append a \t to the end of the data and call the super parse
        rawData = rawData + "\t"
        return super().Parse(headers, rawData)

class PaladinFmtSAMResults(FmtSAMResults):
    def __init__(self, tool, rawData, card=None):
        headers = "ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual	Optionals"
        data = self.Parse(headers, rawData)
        super().__init__(tool, rawData, headers, data)
    
    def Parse(self, headers, rawData):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: ReadName	Flags	RefName	ForwardBase	MapQual	CIGAR	MateRef	MateBase	FragLength	ReadSeq	ReadQual	Optionals
        #returns a dictionary of columnName:data
        #the refname is different, so just change it accordingly
        data = super().Parse(headers, rawData)
        if (data["RefName"] != "*"):
            data["RefName"] = data["RefName"].replace("ARO:", "||ARO:")
        data["ReadName"] = ":".join(data["ReadName"].split(":")[3:])
        return data
    
    def Filter(self, resultList, newReadObj):
        resultList = [r for r in resultList if r.ReadSeq != "*"]
        return super().Filter(resultList, newReadObj)

class KMAResults(Results):
    def __init__(self, tool, rawData, card=None):
        headers = "ReadSeq	NumMapped	Score	Start	Stop	TemplateName	ReadName"
        data = self.Parse(headers, rawData)

        self.ReadSeq = data["ReadSeq"]
        self.read = data["ReadName"]
        self.NumMapped = data["NumMapped"]
        self.Score = data["Score"] 
        self.Start = data["Start"]
        self.Stop = data["Stop"]
        self.TemplateName = data["TemplateName"]
        self.aro = data["TemplateName"].split("|")[4].split(":")[1]
        self.name = data["TemplateName"].split("|")[5].split("[")[0].strip()

        super().__init__(self.read, self.aro, self.name, tool, headers, data)
    
    def Parse(self, headers, rawData):
        #rawdata should be a list or tab sep string
        #headers should be a string
        #headers: ReadSeq NumMapped   Score   Start   Stop TemplateName ReadName
        #returns a dictionary of columnName:data
        metadata = {}
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        col = 0
        for header in headers:
            metadata[header] = data[col].strip()
            col = col + 1
        
        return metadata

    def Filter(self, resultList, newReadObj, evalue = 1e-10, bitscore = 100):
        #resultList is a list of AAFmt6Result Objects
        for i in range(0, len(resultList)):
            if (resultList[i].read == newReadObj.read):
                if (newReadObj.TemplateName != "*"):
                    resultList[i].aro = resultList[i].aro + ";" + newReadObj.aro
                    resultList[i].name = resultList[i].name + ";" + newReadObj.name
        return resultList

class SimTruthResult(Results):
    def __init__(self, tool, rawData, card=None):
        headers = "ReadName\tARO"
        data = self.Parse(headers, rawData, card)

        self.read = data["ReadName"]
        self.aro = data["ARO"]
        self.name = data["Name"]

        super().__init__(self.read, self.aro, self.name, tool, headers, data)
    
    def Parse(self, headers, rawData, card):
        #headers: "ReadName\tARO"
        #returns a dictionary of columnName:data
        card = card.GetCARDALL()
        metadata = {}
        data = rawData.split("\t")
        if not (isinstance(headers, list)):
            headers = headers.split("\t")
        col = 0
        for header in headers:
            metadata[header] = data[col].strip()
            col = col + 1
        if (metadata["ARO"].find(";") > -1):
            #multiple ARO assigned
            names = []
            for aro in metadata["ARO"].split(";"):
                if (aro in card.keys()):
                    names.append(card[aro].name)
                else:
                    names.append("NA")
            metadata["Name"] = ";".join(names)
        else:
            if (metadata["ARO"] in card.keys()):
                metadata["Name"] = card[metadata["ARO"]].name
            else:
                metadata["Name"] = "NA"
        return metadata

class ResultFiles():
    class ResultFileParams():
        def __init__(self, path, methodType, toolName, parameterLabel, classToUse, data) -> None:
            self.path = path
            self.methodType = methodType
            self.toolName = toolName
            self.classToUse = classToUse
            self.parameterLabel = parameterLabel
            self.data = data
            
        def __str__(self) -> str:
            return "{}\t{}\t{}".format(self.methodType, self.toolName, self.parameterLabel)
    
    def __init__(self, listOfFiles, card, truthObjList, outputRoot = "out") -> list:
        fileExtensionClassMap = {"out6DNA":Fmt6Results, "out6AA":AAFmt6Result, "mmseq":MMSeqFmt6Result, "sam":FmtSAMResults, "groot":GrootFmtSAMResults, "kma":KMAResults, "hmm":HMMResults, "paladin":PaladinFmtSAMResults, "truth":SimTruthResult}

        for i in listOfFiles:
            print(i)
            d = "d/d/" + i.split("results/")[1]
            md = d.split("/")
            mt = md[2]
            tool = md[3]
            param = md[4].split("_")[1].split(".")[0]
            if (md[4].find("mmseq") > -1):
                ctu = fileExtensionClassMap["mmseq"]
            elif(md[4].find("blastn") > -1):
                ctu = fileExtensionClassMap["out6DNA"]
            elif(md[4].find("out6") > -1):
                ctu = fileExtensionClassMap["out6AA"]
            elif(md[4].find("groot") > -1):
                ctu = fileExtensionClassMap["groot"]
            elif(md[4].find("paladin") > -1):
                ctu = fileExtensionClassMap["paladin"]
            elif(md[4].find("sam") > -1):
                ctu = fileExtensionClassMap["sam"]
            elif(md[4].find("kma") > -1):
                ctu = fileExtensionClassMap["kma"]                   
            elif(md[4].find("tbl") > -1):
                ctu = fileExtensionClassMap["hmm"]       
            else:
                raise exception("unknown filetype")

            if (i.find("bam") != -1):
                pysam.view("-h", "-o" , i + ".sam", i, catch_stdout=False)
                i = i + ".sam"
            
            if (i.find(".gz") != -1):
                with gzip.open(i, 'r') as f_in, open(i.replace(".gz",""), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                i = i.replace(".gz", "")
                
            #print(datetime.datetime.now())
            #num_lines = sum(1 for line in open(i))
            #print(datetime.datetime.now())

            #count = 0
            #create dictionary of all the reads as keys
            id = []
            with open(i, "r") as fh:
                for line in fh:
                    if (line.startswith('#') == False and line.startswith('@') == False):
                        d = ctu(tool, line, card = card)
                        id.append(d.CleanReadName(d.read))
                   # count = count + 1
                    #if (count % 10000 == 0):
                    #    break
                    #    print(str(count) + "/" + str(num_lines))
            data = dict.fromkeys(set(sorted(id)))

            #print(datetime.datetime.now())
            #count = 0
            with open (i, 'r') as fh:
                for line in fh:
                    if (line.startswith('#') == False and line.startswith('@') == False):
                        d = ctu(tool, line, card = card)
                        readID = d.CleanReadName(d.read)

                        if data[readID] == None:
                            data[readID] = [d]
                        
                        if len(data[readID]) == 0:
                            data[readID] = [d]

                        #non-none value for key at this point, so its safe to just call the class filter function.
                        data[readID] = data[readID][0].Filter(data[readID], d)

                    #count = count + 1
                    #if (count % 10000 == 0):
                        #print(str(count) + "/" + str(num_lines))
                    #    break

            #output it
            outFilePath = outputRoot + "/" + tool + ".tsv"
            with open(outFilePath, "a") as fh:
                for key in data.keys():
                    if (data[key] != None):
                        if (len(data[key]) != 0):
                            result = self.ResultFileParams(i, mt, tool, param, ctu, data[key])
                            toolData = (str(result))
                            for d in result.data: #dump the data
                                d = d.AppendTruthARO(d, truthObjList)
                                d = d.AppendAROFamily(d, card)
                                fh.write(toolData + "\t" + str(d) +  "\n")
    

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--result_file",
        type=str,
        default="../results/files.txt",
        help="text file containing all file paths of result outputs",
    )
    parser.add_argument(
        "-d",
        "--card_dna",
        type=str,
        default="../database/nucleotide_fasta_protein_homolog_model.fasta",
        help="card nucleotide fasta filepath",
    )
    parser.add_argument(
        "-p",
        "--card_protein",
        type=str,
        default="../database/protein_fasta_protein_homolog_model.fasta",
        help="card protein fasta filepath",
    )
    parser.add_argument(
        "-i",
        "--card_index",
        type=str,
        default="../database/aro_index.tsv",
        help="card aro_index.tsv filepath",
    )
    parser.add_argument(
        "-hmmdna",
        "--hmm_dna_cluster_path",
        type=str,
        default="../hmm_nt",
        help="hmm dna cluster filepath",
    )
    parser.add_argument(
        "-hmmprot",
        "--hmm_protein_cluster_path",
        type=str,
        default="../hmm_aa",
        help="hmm dna cluster filepath",
    )
    parser.add_argument(
        "-t",
        "--truth_read_map",
        type=str,
        default="../amr_metagenome.tsv",
        help="filepath to the read > Arg mapping tsv",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./out/",
        help="output path",
    )

    args = parser.parse_args()

    if (os.path.isdir(args.output)):
        shutil.rmtree(args.output)
    os.mkdir(args.output)

    #load CARD
    card = CARD(args.card_dna, args.card_protein, args.card_index)
    card.AssignHMMClusters(args.hmm_dna_cluster_path, args.hmm_protein_cluster_path)  #assign hmm cluster names

    #parse the truth TSV 
    truthDict = {}
    for element in [SimTruthResult("truth",line, card=card) for line in open(args.truth_read_map, 'r')]:
        truthDict[element.read] = element

    #parse the tool results
    ResultFiles([line.strip() for line in open(args.result_file, 'r')], card, truthDict, args.output)
main()