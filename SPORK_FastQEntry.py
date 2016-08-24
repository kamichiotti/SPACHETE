#Fastq entry class

#Imports
import sys

class FastQEntry(object):
    __slots__ = ["read_id","seq","plus_line","quality"]

    def __init__(self,read_id,seq,plus_line,quality):
        """
        Goal: initialize the FastQEntry object with the appropriate info
        Arguments:
            all of the fastq lines one at a time
            read_id, seq, plus_line, and quality

        Returns:
            nothing
        """
        self.read_id = read_id
        self.seq = seq
        self.plus_line = plus_line
        self.quality = quality
        self.clean()

    def get_edge_thirds(self,min_third=20):
        """
        Goal: split this FastQEntry into a 5' and 3' FastQEntry
        Arguments:
            none

        Returns:
            a tuple of FastQEntry objects [5',3']
        """
        third_len = len(self.seq)/3
        if third_len < min_third:
            #sys.stderr.write("Skipping read too short to split in thirds: "+self.read_id+"\n")
            return None,None
        five_prime_seq = self.seq[:third_len]
        three_prime_seq = self.seq[2*third_len:]
        five_prime_read = FastQEntry(self.read_id+"/5_prime",five_prime_seq,self.plus_line,self.quality[:third_len])
        three_prime_read = FastQEntry(self.read_id+"/3_prime",three_prime_seq,self.plus_line,self.quality[2*third_len:])

        return five_prime_read,three_prime_read

    def get_first_last_n(self,third_len=36):
        """
        Goal: very similar to the edge 
        Arguments:
            optional length of n to take (defaulted at 36)

        Returns:
            a tuple of 5' and 3' FastQEntry objects
        """
        #Check to make sure can at least get the first and last third in length
        if len(self.seq) <= third_len*2:
            #sys.stderr.write("Skipping read too short to split in thirds: "+self.read_id+" len = "+str(len(self.seq))+"\n")
            return None,None
        five_prime_seq = self.seq[:third_len]
        three_prime_seq = self.seq[-third_len:]
        five_prime_read = FastQEntry(self.read_id+"/5_prime",five_prime_seq,self.plus_line,self.quality[:third_len])
        three_prime_read = FastQEntry(self.read_id+"/3_prime",three_prime_seq,self.plus_line,self.quality[-third_len:])

        return five_prime_read,three_prime_read

        
    def clean(self):
        """
        Goal: clean up the read id, sequence, plus line, and quality
        Arguments:
            none

        Returns:
            nothing
        """
        #self.read_id = self.read_id.replace(" ","_").replace("\t","_").replace("\n","")
        self.read_id = self.read_id.replace("\n","")
        self.seq = self.seq.replace("U","T").replace("\n","")
        self.plus_line = self.plus_line.replace(" ","_").replace("\n","")
        self.quality = self.quality.replace("\n","")
    
    def __str__(self):
        """
        Goal: return an easy to print string of a FastQEntry
        Arguments:
            none

        Returns:
            An output string
        """
        ret_str = ""
        ret_str += self.read_id+"\n"
        ret_str += self.seq+"\n"
        ret_str += self.plus_line+"\n"
        ret_str += self.quality+"\n"
        return ret_str

    def __lt__(self,other):
        """
        Goal: allow comparison between two FastQEntry objects based on read_id
        Arguments:
            other is a FastQEntry to compare to

        Returns:
            True if the read id of this object is 'less' than that of the other
        """
        return self.read_id < other.read_id


