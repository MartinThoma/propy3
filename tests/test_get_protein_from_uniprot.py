# Core Library
import os
from tempfile import mkstemp

# First party
from propy.GetProteinFromUniprot import GetProteinSequence, GetProteinSequenceFromTxt


def test_main():
    _, result_filepath = mkstemp(suffix="result.txt", prefix="propy3")
    _, target_filepath = mkstemp(suffix="target.txt", prefix="propy3")
    with open(result_filepath, "wb") as savefile:
        with open(target_filepath, "r") as localfile:
            for index, i in enumerate(localfile):
                itrim = i.strip()
                if itrim == "":
                    continue
                else:
                    temp = GetProteinSequence(itrim)
                    print("--------------------------------------------------------")
                    print("The %d protein sequence has been downloaded!" % (index + 1))
                    print(temp)
                    savefile.write(temp + "\n")
                    print("--------------------------------------------------------")

    flag = GetProteinSequenceFromTxt(
        "/home/orient/ProPy/", target_filepath, result_filepath
    )
    print(flag)

    # Cleanup
    os.remove(result_filepath)
    os.remove(target_filepath)
