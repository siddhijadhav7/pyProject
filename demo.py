from threading import Thread
import time
class UIFirewall:
    IPs=[]
    MACs=[]
    def __init__(self):
        obj=Thread(target=self.UI_loop())
        obj.start()

    def UI_loop(self):
        while(True):
            choice=int(raw_input("\n1:Allow_IP\n2:Block_IP\n3:Allow_MAC\n4:Block_MAC\n5:Exit\nEnter your choice:"))
            
            if choice==1:
                IP_address1=raw_input("\nEnter IP address:")
                self.Allow_IP(IP_address1)
            elif choice==2:
                IP_address2=raw_input("\nEnter IP address:")
                self.Block_IP(IP_address2)
            elif choice==3:
                MAC_address1=raw_input("\nEnter MAC address:")
                self.Allow_MAC(MAC_address1)
            elif choice==4:
                MAC_address2=raw_input("\nEnter MAC address:")
                self.Block_MAC(MAC_address2)
            elif choice==5:
                break
            else:
                print "Invalid choice"

    def Allow_IP(self,IP_address1):
        #print '%s'%(IP_address1)
        self.check_IP_add(IP_address1)
        
    def Block_IP(self,IP_address2):
        #print '%s'%(IP_address2)
        self.check_IP_delete(IP_address2)
        
    def Allow_MAC(self,MAC_address1):
        self.check_MAC_add(MAC_address1)
       
    def Block_MAC(self,MAC_address2):
        self.check_MAC_delete(MAC_address2)

    def check_IP_add(self,IP):
        try:
            index=UIFirewall.IPs.index(IP)
            print "%s is already exist"%(IP)
            
        except:
            UIFirewall.IPs.append(IP)

    def check_IP_delete(self,MAC):
        try:
            index=UIFirewall.IPs.index(IP)
            UIFirewall.IPs.remove(IP)
            print "%s IP is blocked now"%(IP)
            
        except:
            print "%s IP is already blocked"%(IP)

    def check_MAC_add(self,MAC):
        try:
            index=UIFirewall.IPs.index(MAC)
            print "%s is already exist"%(MAC)
            
        except:
            UIFirewall.IPs.append(IP)

    def check_IP_delete(self,MAC):
        try:
            index=UIFirewall.IPs.index(MAC)
            UIFirewall.IPs.remove(MAC)
            print "%s IP is blocked now"%(MAC)
            
        except:
            print "%s IP is already blocked"%(MAC) 
            
    def remove_IP(self,IP_address1):
	   try:     
        #print '%s'%(IP_address1)
	    self.check_IP_add(IP_address1)
       except:
        print "%s IP is already deleted"%(IP_address1)    

thread1=UIFirewall()
