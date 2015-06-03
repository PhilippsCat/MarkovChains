class MSM:
    
    def __init__(self,T):
        self.T = T
        self.statdist = self.stationarydistribution()
        print "Open the 'help'-method if you need help"
        
    def help(self):
        from Tkinter import *

        root = Tk()
        frame = Frame(root)
        frame.pack()

        text = Text(root)
        text.insert(INSERT, "This is the documentation of how to use the MSM-class \n next line")
        text.insert(END, "Bye Bye.....")
        text.pack()


root.mainloop()
        
    def stationarydistribution(self):
        EW,EV = linalg.eig(transpose(T))
        statdist = EV[:,0] / sum(EV[:,0])
        return statdist
    
    def eigenvalues(self):
        EW, EV = linalg.eig(T)
        return EW
    
    def timescales(self):
        EW = self.eigenvalues()
        t = []
        for i in range(0,len(EW)):
            t.append(-1 / log(EW[i]))
        return t
    
    def TPT(self, set):
        qminus = np.ones(len(T))
        qplus  = np.ones(len(T))
        n = len(qminus)
        f = np.zeros(shape=(n,n)) # discrete prob. curr.
        pi = self.stationarydistribution()
        effprobcur = np.zeros(shape=(n,n))
        FAB = 0
        for i in range(0,n):
            for j in range(0,n):
                if (i != j):

                    f[i,j] = pi[i]*qminus[i]*T[i,j]*qplus[j]
                    f[j,i] = pi[j]*qminus[j]*T[j,i]*qplus[i]
                    effprobcur[i,j] = max(0,(f[i,j]- f[j,i]))

                    if (i in set):
                        FAB = FAB + effprobcur[i,j]

        KAB = FAB / dot(pi,qminus)
        return f, effprobcur, FAB, KAB