import numpy as np

A, C, G, T = 0, 1, 2, 3
int_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}

# scoring = np.array([[1,-1,-1,-1],
#                     [-1,1,-1,-1],
#                     [-1,-1,1,-1],
#                     [-1,-1,-1,1]])


gapPen = -2 

scoring = np.array([[3,-3,-3,-3],
                    [-3,3,-3,-3],
                    [-3,-3,3,-3],
                    [-3,-3,-3,3]])

class AlignmentFinder(object):
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.path = []
        # self.scoreGrid = None

    def findGlobalAlignment(self):
        # We add 1 to the size (0,0) element or the top left corner
        self.scoreGrid = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self.__computeArray__()
        print (self.scoreGrid)
        return self.__traceback__()

    def __computeArray__(self):
        # Fill the score grid(matching grid)
        for i in range(self.seq1.size+1):
            self.scoreGrid[i,0] = i*gapPen
        for j in range(self.seq2.size+1):
            self.scoreGrid[0,j] = j*gapPen
        # Filling using the score matrix 
        for i in range(1, self.seq1.size+1):
            for j in range(1, self.seq2.size+1):
                self.scoreGrid[i,j] = max( self.scoreGrid[i-1, j-1] + self.__getScore__(i, j),
                                    self.scoreGrid[i-1, j] + gapPen,
                                    self.scoreGrid[i, j-1] + gapPen, 0 )
    def __getScore__(self, i, j):
        # we check it fin its score based on the matching(score) matrix we defined first and the i-1&j-1 
        return scoring[self.seq1[i-1], self.seq2[j-1]]
    
    def __getAlignedPair__(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
        return (n1, n2)

    def __traceback__(self):
        # Translate the score(matching) grid to matched pairs or gaps
        alignedPair= []
        i = self.seq1.size
        j = self.seq2.size
        maxIdx = self.scoreGrid.argmax

        while i >0 and j>0:
            if self.scoreGrid[i-1, j-1] + self.__getScore__(i, j) == self.scoreGrid[i,j]:
                alignedPair.append(self.__getAlignedPair__(i, j))
                i -= 1
                j -= 1
            elif self.scoreGrid[i-1, j] + gapPen == self.scoreGrid[i,j]:
                alignedPair.append(self.__getAlignedPair__(i, 0))
                i -= 1
            else:
                alignedPair.append(self.__getAlignedPair__(0, j))
                j -= 1
        while i > 0:
            alignedPair.append(self.__getAlignedPair__(i, 0))
            i -= 1
        while j > 0:
            alignedPair.append(self.__getAlignedPair__(0, j))
            j -= 1
        alignedPair.reverse()
        return alignedPair  

if __name__ == "__main__":

    seq1 = np.array([T,A,C,G,G,G,C,C,C,G,C,T,A,C], dtype=np.int16)
    seq2 = np.array([T,A,G,C,C,C,T,A,T,C,G,G,T,C], dtype=np.int16)
    aligned = AlignmentFinder(seq1, seq2)
    alignedPairs = aligned.findGlobalAlignment()

    pairShaped = np.array( alignedPairs )
    print( pairShaped[:,0])
    print( pairShaped[:,1])
