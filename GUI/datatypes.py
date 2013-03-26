class InfiniteList(object):
    """
    Creates a special list where removing an element just puts it back at the end of the list
    """
    def __init__(self, values):
        """
        Parameters
        ----------
        values : list

        """
        self.values = values
        
    def pop(self):
        """
        Return the first element, then put the first element back at the end of the list
        """
        val1 = self.values[0]
        self.values.pop(0)
        self.values.append(val1)
        return val1
    
    def prepend(self,item):
        """
        Put the item back to the front of the list
        
        Parameters
        ----------
        item : object
            Thing to push to the beginning of the list
            
        """
        self.values.remove(item)
        self.values.insert(0,item)