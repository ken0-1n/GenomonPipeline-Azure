# Loop through all worker nodes, update hosts file and copy ssh public key to it
# The script make the assumption that the node is called %WORKER+<index> and have
# static IP in sequence order

ADMIN_USERNAME=$1
ADMIN_PASSWORD=$2
WORKER_NAME=$3
WORKER_IP=$4

cd /tmp/torque-5.1.1*

echo 'I update host - '$WORKER_NAME >> /tmp/azuredeploy.log.$$ 2>&1
echo $WORKER_IP $WORKER_NAME >> /etc/hosts
echo $WORKER_IP $WORKER_NAME >> /tmp/hosts.$$
sudo -u $ADMIN_USERNAME sh -c "sshpass -p '$ADMIN_PASSWORD' ssh-copy-id $WORKER_NAME"

# Push packages to compute nodes

sudo -u $ADMIN_USERNAME scp /tmp/hosts.$$ $ADMIN_USERNAME@$WORKER_NAME:/tmp/hosts >> /tmp/azuredeploy.log.$$ 2>&1
sudo -u $ADMIN_USERNAME scp torque-package-mom-linux-x86_64.sh $ADMIN_USERNAME@$WORKER_NAME:/tmp/. >> /tmp/azuredeploy.log.$$ 2>&1
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo -kS sh -c 'cat /tmp/hosts>>/etc/hosts'"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo -kS /tmp/torque-package-mom-linux-x86_64.sh --install"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo -kS /usr/local/sbin/pbs_mom"
sudo sh -c "echo '$WORKER_NAME' >> /var/spool/torque/server_priv/nodes"


