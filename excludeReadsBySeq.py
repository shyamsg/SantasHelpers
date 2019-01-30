import sys


if (len(sys.argv) != 4):
    print "Usage: excludeReadsBySeq.py NegativeControlFastq TargetFastq OutputFastq"
    sys.exit(1)

negControlFile = open(sys.argv[1])
seqsToExclude = []
cnt = 0
for line in negControlFile:
    cnt = cnt + 1
    if cnt%4 == 2:
        seqsToExclude.append(line.strip())
negControlFile.close()
print "Loaded", len(seqsToExclude), "sequences from negative control fastq file."

TargetFastq = open(sys.argv[2])
OutputFastq = open(sys.argv[3], "w")
quartet = ["","","",""]
cnt = 0
outcnt = 0
for line in TargetFastq:
    quartet[cnt%4] = line.strip()
    if cnt%4 == 3:
        if quartet[1] not in seqsToExclude:
            OutputFastq.write("\n".join(quartet)+"\n")
            outcnt = outcnt + 1
    cnt = cnt + 1

TargetFastq.close()
OutputFastq.close()
print "Processed", int(cnt/4), "reads from input fastq, excluded", int(cnt/4) - outcnt, "reads, retained", outcnt, "reads."
