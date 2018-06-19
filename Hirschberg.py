import numpy as np

A, C, G, T = 0, 1, 2, 3
int_to_char = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

gapPen = -1
# scoring = np.array([[10, -3, -1, -4],
#                     [-3, 9, -5, 0],
#                     [-1, -5, 7, -3],
#                     [-4, 0, -3, 8]])

scoring = np.array([[3,-3,-3,-3],
                    [-3,3,-3,-3],
                    [-3,-3,3,-3],
                    [-3,-3,-3,3]])


def scoreCalc( seq1, seq2 ):
    # This calculate the of the scoreGrid for those sequences seq1 and seq2
    scoreGrid = np.zeros((seq1.size+1, seq2.size+1), dtype=np.int16)
    # Fill the score grid(matching grid)
    for i in range(seq1.size+1):
        scoreGrid[i, 0] = i*gapPen
    for j in range(seq2.size+1):
        scoreGrid[0, j] = j*gapPen
    # Filling using the score matrix
    for i in range(1, seq1.size+1):
        for j in range(1, seq2.size+1):
            scoreGrid[i, j] = max(scoreGrid[i-1, j-1] + getScore(seq1, seq2, i, j),
                                  scoreGrid[i-1, j] + gapPen,
                                  scoreGrid[i, j-1] + gapPen)
    
    return scoreGrid

def NWscore(seq1, seq2):
    # This calculate the last row of the of the scoreGrid for those sequences seq1 and seq2
    scoreGrid = scoreCalc( seq1, seq2 )
    return scoreGrid[seq1.size, :]


def getScore(seq1, seq2, x, y):
    # we check it fin its score based on the matching(score) matrix we defined first and the i-1&j-1
    return scoring[seq1[x-1], seq2[y-1]]


def getAlignedPair(seq1, seq2, i, j):
    n1 = int_to_char[seq1[i-1]] if i > 0 else '_'
    n2 = int_to_char[seq2[j-1]] if j > 0 else '_'
    return (n1, n2)


def traceback(seqX, seqY):
    # Translate the score(matching) grid to matched pairs or gaps
    alignedPair = []
    i = seqX.size
    j = seqY.size
    # scoreGrid = scoreCalc( seqX, seqY )
    if j == 0:
        # This is search for direct mismatch and add gaps for the second sequence
        while i > 0:
            alignedPair.append(getAlignedPair(seqX, seqY, i, 0))
            i -= 1
    elif i == 0:
        # This is search for direct mismatch and add gaps for the first sequence
        while j > 0:
            alignedPair.append( getAlignedPair(seqX, seqY,0, j))
            j -= 1
    elif i == 1 or j == 1:
        # This is just normal Needleman-Wunsch Algorithm 
        scoreGrid = scoreCalc( seqX, seqY )
        while i > 0 or j > 0:
            if i > 0 and j > 0 and scoreGrid[i-1, j-1] + getScore( seqX,seqY, i, j) == scoreGrid[i, j]:
                alignedPair.append(getAlignedPair( seqX,seqY, i, j))
                i -= 1
                j -= 1
            elif i > 0 and scoreGrid[i-1, j] + gapPen == scoreGrid[i, j]:
                alignedPair.append( getAlignedPair( seqX,seqY, i, 0))
                i -= 1
            else:
                alignedPair.append(getAlignedPair(seqX, seqY, 0, j))
                j -= 1
    else:
        # Here is were all the magic happens, and by magic I mean Divide&Conquer happens
        # Where the First sequence is divided
        iMid = int(i/2)
        scoreL = NWscore( seqX[0:iMid], seqY )
        scoreR = NWscore( seqX[iMid:i][::-1], seqY[::-1] )
        sum = np.add(scoreL, scoreR[::-1] )
        
        jMid = np.argmax( sum )
        alignedPair+= traceback ( seqX[0:iMid], seqY[ 0: jMid ] )
        alignedPair+= traceback ( seqX[iMid:i], seqY[ jMid: j ] )  
    
    print(alignedPair)
    return alignedPair


if __name__ == "__main__":
    # array1 = ['G', 'T', 'A', 'C', 'A', 'G', 'T','A']
    # array2 = ['G', 'G', 'T', 'A', 'C', 'G', 'T']
    seq1 = np.array([A,G,T,A,C,G,C,A], dtype=np.int16)
    seq2 = np.array([T,A,T,G,C], dtype=np.int16)
    aligned = traceback( seq1, seq2 )
    # print(aligned)
    # aligned.reverse()
    # print(aligned)
    pairShaped = np.array(aligned)
    print(pairShaped)
    # print(pairShaped[:, 1])
